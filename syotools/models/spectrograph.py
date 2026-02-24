#!/usr/bin/env python
"""
Created on Sat Oct 15 16:56:40 2016

@author: gkanarek, tumlinson
"""

import numpy as np
import astropy.units as u
from astropy.table import QTable

from syotools.models.base import PersistentModel
from syotools.models.source_exposure import SourceSpectrographicExposure
from syotools.spectra.utils import mirror_efficiency, set_coating
from syotools.defaults import default_spectrograph, default_spectropolarimeter
from hwo_sci_eng.utils import read_yaml 

class Spectrograph(PersistentModel):
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

    def __init__(self, default_model = default_spectrograph, **kw):
        self.telescope = None
        self.exposures = []

        self._lumos_default_file = ''

        self.name = ''
        self.modes = []
        self.descriptions = {}
        self.bef = np.zeros(0, dtype=float) * (u.erg / u.cm**2 / u.s / u.pix)
        self.R = 0. * u.dimensionless_unscaled
        self.wave = np.zeros(0, dtype=float) * u.AA
        self.aeff = np.zeros(0, dtype=float) * u.cm**2
        self.wrange = np.zeros(2, dtype=float) * u.AA
        self._mode = ''
        super().__init__(default_model, **kw)


    #Property wrapper for mode, so that we can use a custom setter to propagate
    #mode updates to all the rest of the parameters

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, new_mode):
        """
        Mode is used to set all the other parameters
        """ 

        nmode = new_mode.upper()
        if self._mode == nmode or nmode not in self.modes:
            return
        self._mode = nmode
        table = QTable.read(self._lumos_default_file, nmode)
                
        self.R = table.meta['R'] * u.pix
        self.wave = table['Wavelength'].copy()  # Copy to remove FITS file weakrefs
        self.bef = table['BEF'].copy()          # Copy to remove FITS file weakrefs
        self.aeff = table['A_Eff'].copy()       # Copy to remove FITS file weakrefs
        wrange = np.array([table.meta['WAVE_LO'], table.meta['WAVE_HI']]) * u.AA
        self.wrange = wrange

    @property
    def delta_lambda(self):
        wave, R = self.recover('wave', 'R')
        return wave / R
    
    def create_exposure(self, source=None):
        new_exposure = SourceSpectrographicExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure

    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.spectrograph = self
        exposure.telescope = self.telescope
        exposure.calculate()

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

