#!/usr/bin/env python
"""
Created on Fri Oct 14 21:31:18 2016
@author: gkanarek, tumlinson
"""
import copy

import numpy as np
import astropy.constants as const
import astropy.units as u
import synphot as syn
from synphot.models import Empirical1D

from syotools.models.base import PersistentModel
from syotools.spectra.utils import mag_from_sed, mirror_efficiency, set_coating
from syotools.utils.yaml_utils import simplify_data, complexify_data
from hwome.core.navigator import DataModel


class Instrument(PersistentModel):
    """
    The basic instrument class, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with the camera
        exposures    - the list of Exposures taken with this camera

        _default_model - used by PersistentModel

    The following are attributes I haven't included, and the justification:
        R_effective - this doesn't seem to be used anywhere
    """

    def __init__(self, telescope, **kw):

        self.telescope = telescope
        self.exposures = []
        self.name = ''
        #self.pivotwave = np.zeros(1, dtype=float) * u.nm
        self.bandnames = ['']
        self.channels = [([],0)]
        self.fiducials = np.zeros(1, dtype=float) * u.nm
        self.total_qe = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.ap_corr = np.ones(1, dtype=float) * u.dimensionless_unscaled
        self.bandpass_r = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.dark_current = np.zeros(1, dtype=float) * (u.electron / u.s / u.pixel)
        #self.detector_rn = np.zeros(1, dtype=float) * (u.electron / u.pixel)**0.5
        self.sky = syn.spectrum.SourceSpectrum(Empirical1D, points=[0.1,10000, 20000] << u.AA, lookup_table=[22,22,22] << u.ABmag) # Hardcode a 22nd magnitude background
        #super().__init__(default_camera, **kw)

    def diffraction_limit(self, wavelength: u.Quantity) -> u.Quantity:
        """
        Calculate the diffraction limit for a given wavelength.
        """
        effective_diameter = self.recover('telescope.effective_diameter')
        ap_nm = effective_diameter.to(u.nm)
        diff_limit_radians = 1.22 * u.rad * wavelength.to(u.nm) / ap_nm
        return diff_limit_radians.to(u.arcsec)


    @property
    def diff_limit_fwhm(self):
        """
        Diffraction-limited PSF FWHM.
        """

        configuration, effective_diameter = self.recover('configuration',
                                                       'telescope.effective_diameter')
        diff_limit_wavelength = configuration["diffraction_limit"]

        #result = (1.22 * u.rad * diff_limit_wavelength / aperture).to(u.arcsec)
        result = (1.03 * u.rad * diff_limit_wavelength / effective_diameter).to(u.arcsec)
        return result

    # UNFINISHED
    def ee(self,wavelength):
        effective_aperture = self.recover('telescope.effective_aperture')
        a = effective_radius
        k = 2 * np.pi / wavelength
        q = np.arange()
        x = 1 - [J0(np.pi * r)]**2 - [J1(np.pi * r)]**2

    def fwhm_psf(self, wave):
        """
        Calculate the FWHM of the camera's PSF.
        """
        #Convert to Quantity for calculations.
        effective_aperture = self.recover('telescope.effective_diameter')
        configuration, diff_fwhm = self.recover('configuration',
                                             'diff_limit_fwhm')

        diff_limit = configuration["diffraction_limit"]

        #fwhm = (1.22 * u.rad * wave / aperture).to(u.arcsec)
        fwhm = (1.03 * u.rad * wave / effective_aperture).to(u.arcsec)
        
        #only use these values where the wavelength is greater than the diffraction limit
        fwhm = np.where(wave > diff_limit, fwhm.value, diff_fwhm.value) * u.arcsec

        return fwhm


    def _c_thermal(self, wave, verbose=False):
        """
        Calculate the thermal emission counts for the telescope.
        """

        #Convert to Quantities for calculation.
        (diameter, ota_emissivity, configuration) = self.recover(
                'telescope.effective_diameter',  'telescope.ota_emissivity', 'configuration')
        total_qe = configuration["detector"]["total_qe"]
        pixel_scale = configuration["pixel_scale"]


        box = self._sn_box(wave, verbose)

        h = const.h.to(u.erg * u.s) # Planck's constant erg s
        c = const.c.to(u.cm / u.s) # speed of light [cm / s]

        energy_per_photon = h * c / wave.to(u.cm) / u.ph

        D = diameter.to(u.cm) # telescope diameter in cm

        

        pephot = self._planck(wave) / energy_per_photon

        if verbose:
            print('Planck spectrum: {}'.format(self.nice_print(self.planck(wave))))
            print('Planck / E_phot: {}'.format(self.nice_print(pephot)))
            print('E_phot: {}'.format(self.nice_print(energy_per_photon)))
            #print('Omega: {}'.format(self.nice_print(Omega)))

        thermal = (ota_emissivity[0] * self._planck(wave) / energy_per_photon *
    			(np.pi / 4. * D**2 * u.AA**-1))

        # omega is the size of the extraction box in steradians
        Omega = (pixel_scale**2 * self._sn_box(wave, False)).to(u.sr)
        thermal *= Omega

        thermal = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=thermal.value * syn.units.PHOTLAM)

        return thermal

    def _planck(self, wave):
        """
        Planck spectrum for the various wave bands.
        """
        #Convert to Quantities for calculation
        configuration = self.recover('configuration')
        temperature = configuration["detector"]["thermal"]

        wave = wave.to('cm')
        if isinstance(temperature, u.Quantity):
            temp = temperature
        else:
            temps = temperature[0] * u.Unit(temperature[1])
            temp = temps.to('K')
        h = const.h.to(u.erg * u.s) # Planck's constant erg s 
        c = const.c.to(u.cm / u.s) # speed of light [cm / s] 
        k = const.k_B.to(u.erg / u.K) # Boltzmann's constant [erg deg K^-1] 
        x = 2. * h * c**2 / wave**5 
        exponent = (h * c / (wave * k * temp)) 
    
        result = (x / (np.exp(exponent)-1.)).to(u.erg / u.s / u.cm**3) / u.sr

        return result


    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.instrument = self
        exposure.telescope = self.telescope
        exposure.calculate()

    def set_from_hwome(self, channelname, ins_type):
        self.configuration = {}
        try:
            instrument, channel = channelname.split(".")
            instrument_data = getattr(self.telescope.hwo_data, instrument)
        except ValueError:
            raise ValueError("Need a name + channel (e.g. 'HRI-S.NIR')")
        except KeyError:
            raise KeyError(f"Unrecognized Instrument {instrument}.\n Legal values are {self.telescope.hwo_data.Instrument.name}")

        instrument_data = self.telescope.hwo_data[instrument]

        self.name = channel
        instrument_data.Channel
        try:
            channel_data = getattr(instrument_data, channel)
        except KeyError:
            raise KeyError(f"Unrecognized Channel {channel}.\n Legal values are {instrument_data.Channel.name.keys()}")
        channel_data = instrument_data[channel]


        try:
            channel_filters = list(channel_data.Filter.name.keys())
        except TypeError:
            channel_filters = [channel_data.Filter.name.value]
        self.configuration["band"] = {}
        self.configuration["channel_filters"] = []

        for channel_filter in channel_filters:
            filter_name = channel_filter.split(".")[-1]
            fancy_name = ".".join((channel, filter_name))

            # only load optical elements if they're the correct type
            try:
                getattr(channel_data[filter_name], "Grating")
                kind = "disperser"
                if "ifu" in filter_name.lower():
                    kind = "ifs"
            except KeyError:
                kind = "filter"
            if ins_type.lower() in ("imager") and kind != "filter":
                continue
            elif ins_type.lower() in ("spectrograph") and kind != "disperser":
                continue
            elif ins_type.lower() in ("ifs") and kind != "ifs":
                continue

            self.configuration["channel_filters"].append(fancy_name)
            # this commands HWOME to walk down the entire optical path of the telescope down to the filter(grating) and collect all of the optics.
            thru = self.telescope.hwo_data.OpticalPath.select(instrument=instrument, channel=channel, filter = filter_name).throughput(include_detector=False)
            # then we multiply all of them together
            total_throughput = np.prod(thru.q, axis=0)

            # and store for later retrieval
            band = self.load_throughput(thru.w, total_throughput)
            wavemin = band.avgwave() - band.rectwidth()/2
            wavemax = band.avgwave() + band.rectwidth()/2
            self.configuration["band"][fancy_name] = {"internal_name": filter_name, "bandpass": band, "original_wave": thru.w, "original_thru": total_throughput, "effective_wavelength": band.avgwave(), 
                                                    "wave_min": wavemin, "wave_max": wavemax, "optics": len(thru.value.keys())}
            if kind in ("disperser", "ifs"):
                grating_resolution = channel_data[filter_name].Grating.spectral_resolution.q
                self.configuration["band"][fancy_name]["resolution"] = float(grating_resolution)
            self.configuration["band"][fancy_name]["kind"] = kind

        self.configuration["diffraction_limit"] = channel_data.diffraction_limited.q
        self.configuration["pixel_scale"] = channel_data.plate_scale.q
        self.configuration["fov_x"] = channel_data.fov_x.q
        self.configuration["fov_y"] = channel_data.fov_y.q

        self.configuration["detector"] = {}
        for detector in channel_data.Detector:
            self.configuration["detector"]["internal_name"] = detector.name
            self.configuration["detector"]["read_noise"] = detector.read_noise.q
            self.configuration["detector"]["thermal"] = detector.temperature.q
            self.configuration["detector"]["dark_current"] = detector.dark_current.q / u.pix #* u.electron / u.pix**2 / u.ct # needs to be electrons per pixel per second
            w = detector.qe.w
            t = detector.qe.q
            self.configuration["detector"]["total_qe"] = self.load_throughput(w, t)
            self.configuration["detector"]["original_qe_wave"] = w
            self.configuration["detector"]["original_qe_thru"] = t

    def load_throughput(self, wave, thru):
        return syn.spectrum.SpectralElement(Empirical1D, points=wave, lookup_table=thru)

    def load_from_dict(self, config):
        """
        Restore an instrument from a stored dictionary
        """
        config = complexify_data(config)

        for band in config["band"]:
            config["band"][band]["bandpass"] = self.load_throughput(config["band"][band]["original_wave"], config["band"][band]["original_thru"])

        config["detector"]["total_qe"] = self.load_throughput(config["detector"]["original_qe_wave"], config["detector"]["original_qe_thru"])

        self.configuration = config

    def save_to_dict(self):
        """
        Save an instrument to a dictionary - scrub the Synphot objects so we don't have to flatten them

        Returns
        -------
        config: dict
            A configuration dictionary
        """
        config = copy.deepcopy(self.configuration)
        for band in config["band"]:
            del config["band"][band]["bandpass"]
        
        del config["detector"]["total_qe"]

        config = simplify_data(config)

        return config

    def set_from_sei(self, name): 
        """
        This method is implemented at the subclass level, not here
        (so in Camera, Spectrograph, IFS, etc.)
        """
        raise NotImplementedError
