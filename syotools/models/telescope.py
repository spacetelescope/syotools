#!/usr/bin/env python
"""
Created on Fri Oct 14 20:28:51 2016

@author: gkanarek, tumlinson
"""

from syotools.models.base import PersistentModel
from syotools.defaults import default_telescope
from syotools.utils import pre_encode
from syotools.utils.jsonunit import str_jsunit
import astropy.units as u #for unit conversions
import numpy as np
from syotools.sci_eng_interface import read_json 
from hwo_sci_eng.utils import read_yaml 

# would be true with freestanding Sci-Eng-Interface 

class Telescope(PersistentModel):
    """
    The basic telescope class, which provides parameter storage for 
    optimization.
    
    Attributes: #adapted from the original in Telescope.py
        name - The name of the telescope (string)
        aperture - The size of the primary telescope aperture, in meters (float)
        unobscured_fraction - The fraction of the primary mirror which is not obscured (float)
        temperature - instrument temperature, in Kelvin (float)
        ota_emissivity - emissivity factor for a TMA (float)
        diff_limit_wavelength - diffraction limit wavelength, in nm (float)
        
        _default_model - used by PersistentModel
        
        cameras - the Camera objects for this telescope
    """
    _default_model = default_telescope
    
    cameras = []
    spectrographs = []
    coronagraphs = [] 
    
    name = ''
    aperture = pre_encode(0. * u.m)
    temperature = pre_encode(0. * u.K)
    ota_emissivity = pre_encode(0. * u.dimensionless_unscaled)
    diff_limit_wavelength = pre_encode(0. * u.nm)
    unobscured_fraction = pre_encode(1. * u.dimensionless_unscaled)

    verbose = False 
        
    @property
    def diff_limit_fwhm(self):
        """
        Diffraction-limited PSF FWHM.
        """
        
        diff_limit_wavelength, aperture = self.recover('diff_limit_wavelength',
                                                       'aperture')
        
        #result = (1.22 * u.rad * diff_limit_wavelength / aperture).to(u.arcsec)
        result = (1.03 * u.rad * diff_limit_wavelength / aperture).to(u.arcsec)
        return pre_encode(result)
    
    @property
    def effective_aperture(self):
        unobscured, aper = self.recover('unobscured_fraction', 'aperture')
        return pre_encode(np.sqrt(unobscured) * aper)
    
    def add_camera(self, camera):
        self.cameras.append(camera)
        camera.telescope = self
    
    def add_spectrograph(self, spectrograph):
        self.spectrographs.append(spectrograph)
        spectrograph.telescope = self

    def add_coronagraph(self, coronagraph):
        self.coronagraphs.append(coronagraph)
        coronagraph.telescope = self

    def set_from_json(self,name): 
        if self.verbose: 
            print('Setting Telescope to: ', name) 
        
        if ('EAC2' in name): tel = read_json.eac2()
        if ('EAC3' in name): tel = read_json.eac3()
        
        self.name = tel['name'] 
        self.aperture = tel['aperture_od'] * u.m 
        self.temperature = tel['temperature_K'] * u.K 
        self.diff_limited_wavelength = tel['diff_limited_wavelength'] * u.nm 
        self.unobscured_fraction = tel['unobscured_fraction'] 

    def set_from_yaml(self,name): 

        if ('EAC1' in name): tel = read_yaml.eac1()

        self.name = tel['name'] 
        self.inscribing_aperture = tel["PM_aperture"]['segmentation_parameters']['inscribing_diameter'][0] * u.m 
        self.circumscribing_aperture = tel["PM_aperture"]['segmentation_parameters']['circumscribing_diameter'][0] * u.m 

        if ('EAC2' in name): print("You can't set EAC2 from the SEI YAML file") 
        if ('EAC3' in name): print("You can't set EAC3 from the SEI YAML file") 

       



