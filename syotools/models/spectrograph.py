#!/usr/bin/env python
"""
Created on Sat Oct 15 16:56:40 2016

@author: gkanarek, tumlinson
"""

import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.table import QTable
import synphot as syn
import stsynphot as stsyn
from synphot.models import Empirical1D

from .instrument import Instrument
from syotools.models.source_exposure import SourceSpectrographicExposure
from syotools.spectra.utils import mirror_efficiency, set_coating
from syotools.defaults import default_spectrograph, default_spectropolarimeter
from hwome.core.navigator import DataModel

class Spectrograph(Instrument):
    """
    The basic spectrograph class, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with this spectrograph
        exposures    - the list of Exposures taken with this spectrograph

        name         - name of the spectrograph (string)

        modes        - supported observing modes (list)
        descriptions - description of supported observing modes (dict)
        mode         - current observing mode (string)
        bef          - background emission function in erg/cm2/pixel/s ement (float array)
        R            - spectral resolution (float)
        wrange        - effective wavelength range (2-element float array)
        wave         - wavelength in Angstroms (float array)
        aeff         - effective area at given wavelengths in cm^2 (float array)

        _lumos_default_file - file path to the fits file containing LUMOS values

        _default_model - used by PersistentModel
    """

    def __init__(self, telescope, default_model = default_spectrograph, **kw):
        self.telescope = telescope
        self.exposures = []

        self._lumos_default_file = ''

        self.name = ''
        self.descriptions = {}
        self.sky = syn.spectrum.SourceSpectrum(Empirical1D, points=[0.1,20000] << u.AA, lookup_table=[24,24] << u.ABmag)
        self.sky = self.sky.normalize(24 * u.ABmag, stsyn.spectrum.band("johnson,v"))
        self.R = 0. * u.dimensionless_unscaled
        self.wave = np.zeros(0, dtype=float) * u.AA
        self.aeff = np.zeros(0, dtype=float) * u.cm**2
        self.wrange = np.zeros(2, dtype=float) * u.AA
        self._band = None
        #super().__init__(default_model, **kw)


    #Property wrapper for band, so that we can use a custom setter to propagate
    #band updates to all the rest of the parameters

    @property
    def n_bands(self):
        return len(self.bands)

    @property
    def band(self):
        return self._band

    @property
    def bands(self):
        return [x for x in self.configuration["band"] if self.configuration["band"][x]["kind"] == "disperser"]

    @band.setter
    def band(self, new_band):
        """
        Band is used to set all the other parameters
        """ 

        nband = new_band.upper()
        if self._band == nband or nband not in self.configuration["channel_filters"]:
            return
        self._band = nband

        self.R = self.configuration["band"][nband]["resolution"]
        self.wave = self.configuration["band"][nband]["bandpass"].waveset
        self.sky = syn.spectrum.SourceSpectrum(Empirical1D, points=self.wave, lookup_table=np.ones_like(self.wave.value) * 24 << u.ABmag)
        self.sky = self.sky.normalize(24 * u.ABmag, stsyn.spectrum.band("johnson,v"))
        self.aeff = self.configuration["band"][nband]["bandpass"]
        wrange = np.array((np.min(self.wave.value), np.max(self.wave.value)))
        self.wrange = wrange

    @property
    def delta_lambda(self):
        wave, R = self.recover('wave', 'R')
        R = R << u.pix # HWOME's definition is unitless
        return wave / R

    def _sn_box(self, wave, verbose):
        """
        Calculate the number of pixels in the SNR computation box.
        """

        Phi = self.configuration["pixel_scale"]
        sn_box = np.round(3. * self.fwhm_psf(wave) / Phi)

        if verbose:
            print('PSF width: {}'.format(self.nice_print(self.fwhm_psf(wave))))
            print('SN box height: {}'.format(self.nice_print(sn_box)))

        return sn_box * 1 * u.pix

    def _print_initcon(self, verbose):
        if verbose: #These are our initial conditions
            print('Telescope diffraction limit: {}'.format(self.telescope.diff_limit_wavelength))
            print('Telescope effective_diameter: {}'.format(self.telescope.effective_diameter))
            print('Instrument temperature: {}'.format(self.configuration["detector"]["thermal"]))
            print("R", self.R)
            print('Resolution: {}'.format(self.nice_print(self.R)))
            print('Pixel sizes: {}'.format(self.configuration["pixel_scale"]))
            #print('AB mag zero points: {}'.format(self.nice_print(self.ab_zeropoint)))
            print('Quantum efficiency: {}'.format(self.nice_print(self.configuration["detector"]["total_qe"].waveset)))
            #print('Aperture correction: {}'.format(self.ap_corr))
            print('Detector read noise: {}'.format(self.configuration["detector"]["read_noise"]))
            print('Dark rate: {}'.format(self.configuration["detector"]["dark_current"]))


    def create_exposure(self, source=None):
        new_exposure = SourceSpectrographicExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure

    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.instrument = self
        exposure.telescope = self.telescope
        exposure.calculate()

    def transform_flux(self, spectrum, wave):
        effective_area = self.recover("telescope.effective_area")
        flux = syn.units.convert_flux(wave, spectrum(wave), u.erg / u.s / u.cm**2 / u.AA)
        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.ct
        return flux / phot_energy * effective_area

    def set_from_sei(self, name):

        if ('uvi' in name.lower()): uvi = read_yaml.uvi()
        
        # the "uvi" dictionary returned by read_yaml is nested, and therefore awkward  
        # when summoning individual entries. And often, we do not need the individual 
        # components. So, we are going to break this dictionary up and carry the 
        # pieces separately:  

        self.fuv_imager = uvi['FUV_Imager']  

        self.fuv_mos = uvi['FUV_MOS'] 
        
        self.nuv_mos = uvi['NUV_MOS']

        self.msa = uvi['MSA'] 

        self.mcp = uvi['MCP']

        self.cmos = uvi['CMOS']

        self.imager_mirrors = {}
        for mirror in range(1, self.fuv_imager["N_refl_optics"][0] + 1):
            self.imager_mirrors[f"mirror{mirror}"] = set_coating(self.fuv_imager)
        self.instrument_efficiency_imager = mirror_efficiency(self.imager_mirrors)

        for disperser in self.fuv_mos:
            fuvmos_mirrors = {}
            if isinstance(self.fuv_mos[disperser], dict):
                for mirror in range(1, self.fuv_mos[disperser]["N_refl_optics"][0] + 1):
                    fuvmos_mirrors[f"mirror{mirror}"] = set_coating(self.fuv_mos[disperser])
            setattr(self, f"instrument_efficiency_{disperser}", mirror_efficiency(fuvmos_mirrors))

        for disperser in self.nuv_mos:
            nuvmos_mirrors = {}
            if isinstance(self.nuv_mos[disperser], dict):
                for mirror in range(1, self.nuv_mos[disperser]["N_refl_optics"][0] + 1):
                    nuvmos_mirrors[f"mirror{mirror}"] = set_coating(self.nuv_mos[disperser])
            setattr(self, f"instrument_efficiency_{disperser}", mirror_efficiency(nuvmos_mirrors))

        # Get a handy list of available dispersers
        self.modes = list(self.fuv_mos.keys())
        self.modes.extend(list(self.nuv_mos.keys()))
        self.modes.remove("G165LL") # No data to support mode
        self.modes.remove("G700L") # No data to support mode

class Spectropolarimeter(Spectrograph):
    """
    The basic spectropolarimeter class for POLLUX, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with this spectrograph
        exposures    - the list of Exposures taken with this spectrograph

        name         - name of the spectrograph (string)

        modes        - supported observing modes (list)
        descriptions - description of supported observing modes (dict)
        mode         - current observing mode (string)
        bef          - background emission function in erg/cm2/s/pixel (float array)
        R            - spectral resolution (float)
        wrange        - effective wavelength range (2-element float array)
        wave         - wavelength in Angstroms (float array)
        aeff         - effective area at given wavelengths in cm^2 (float array)

        _lumos_default_file - file path to the fits file containing LUMOS values

        _default_model - used by PersistentModel
    """

    def __init__(self, **kw):
        super().__init__(default_spectropolarimeter, **kw)

