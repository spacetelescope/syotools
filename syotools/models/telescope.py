#!/usr/bin/env python
"""
Created on Fri Oct 14 20:28:51 2016
@authors: gkanarek, tumlinson
"""

from syotools.models.base import PersistentModel
from syotools.defaults import default_telescope
from syotools.utils import pre_encode
from syotools.utils.jsonunit import str_jsunit
import astropy.units as u 
import numpy as np
import os, yaml
from syotools.sci_eng_interface import read_json 
from hwo_sci_eng.utils import read_yaml 

class Telescope(PersistentModel):
    """
    The basic telescope class to provide telescope parameter storage. 
    
    Attributes: #adapted from the original in Telescope.py
        name - The name of the telescope (string)
        effective_aperture - The size of the primary telescope aperture, in meters (float)
            note: there is no such thing as "aperture", there is only "effective aperture" 
                for a circular/keystone primary, this is just the diameter of the circle. 
                for a hex-pattern segmented primary, this is the diameter of a circle 
                with the same area as the summed area of all the hex segments. 
                all code should use ONLY effective_aperture
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
    effective_aperture = pre_encode(0. * u.m) # held over from before SEI integration 
    temperature = pre_encode(0. * u.K) # held over from before SEI integration
    ota_emissivity = pre_encode(0. * u.dimensionless_unscaled)  # held over from before SEI integration
    diff_limit_wavelength = pre_encode(0. * u.nm)  # held over from before SEI integration
    unobscured_fraction = pre_encode(1. * u.dimensionless_unscaled)  # held over from before SEI integration

    verbose = False 
        
    @property
    def diff_limit_fwhm(self):
        """
        Diffraction-limited PSF FWHM.
        """
        
        diff_limit_wavelength, effective_aperture = self.recover('diff_limit_wavelength',
                                                       'effective_aperture')
        
        #result = (1.22 * u.rad * diff_limit_wavelength / effective_aperture).to(u.arcsec)
        result = (1.03 * u.rad * diff_limit_wavelength / effective_aperture).to(u.arcsec)
        return pre_encode(result)
    
#    @property
#    def effective_aperture(self):
#        self.effective_aperture = 2. * (self.total_collecting_area / np.pi )**0.5
#        return 

#    @effective_aperture.setter
#    def effective_aperture(self): 
#        self.effective_aperture = 2. * (self.total_collecting_area / np.pi )**0.5
#        return
    
    def add_camera(self, camera):
        self.cameras.append(camera)
        camera.telescope = self
    
    def add_spectrograph(self, spectrograph):
        self.spectrographs.append(spectrograph)
        spectrograph.telescope = self

    def add_coronagraph(self, coronagraph):
        self.coronagraphs.append(coronagraph)
        coronagraph.telescope = self

    def hexagon_area(self, side): 
        return 3. * 3.**0.5 / 2. * side**2
    
    def set_coating(self, mirror, coating): 
        """ Sets the reflectivity curve for a mirror surface by 
            adding / modifying that mirror's component dictionary 
            doing this with file I/O here is a bit inelegant but 
            works for the first iteration of this capability. 
            JT 102524
        """
        with open(os.getenv('SCI_ENG_DIR') + '/obs_config/Tel/'+coating+'_refl.yaml', 'r') as f:
            coating_dict = yaml.load(f, Loader=yaml.SafeLoader)

        mirror['coating_name'] = coating 
        mirror['coating_wave'] = coating_dict['wavelength'] * u.nm 
        mirror['coating_refl'] = coating_dict['reflectivity'] * u.dimensionless_unscaled

    def set_from_sei(self, name): 
        if name == 'EAC1': 
            self.set_from_yaml(name)
        elif name == 'EAC2': 
            self.set_from_yaml(name)
        elif name == 'EAC3': 
            self.set_from_yaml(name)
        else: 
            print('We do not have SEI information for: ', name)
            raise NotImplementedError

    def set_from_json(self,name): 
        if self.verbose: 
            print('Setting Telescope to: ', name) 
        
        if ('EAC2' in name): tel = read_json.eac2()
        if ('EAC3' in name): tel = read_json.eac3()
        
        self.name = tel['name'] 
        self.effective_aperture = tel['aperture_od'] * u.m 
        self.temperature = tel['temperature_K'] * u.K 
        self.diff_limited_wavelength = tel['diff_limited_wavelength'] * u.nm 
        self.unobscured_fraction = tel['unobscured_fraction'] 

    def set_from_yaml(self, name): 

        if ('EAC1' in name): tel = read_yaml.eac1()

        if ('EAC2' in name): 
            print("You can't set EAC2 from the SEI YAML file") 
            raise NotImplementedError

        if ('EAC3' in name): 
            print("You can't set EAC3 from the SEI YAML file") 
            raise NotImplementedError
        
        # the "tel" dictionary returned by read_yaml is nested, and therefore awkward  
        # when summoning individual entries. And often, we do not need the individual 
        # mirrors. So, we are going to break this dictionary up and carry the 
        # mirrors and other pieces separately:  

        self.pm = tel['PM_aperture'] #primary 
        self.set_coating(self.pm, 'XeLiF')
        self.sm = tel['SM'] # secondary  
        self.set_coating(self.pm, 'XeLiF')
        self.m3 = tel['M3'] # tertiary 
        self.set_coating(self.pm, 'XeLiF')
        self.m4 = tel['M4'] # fold mirror (?) 
        self.set_coating(self.pm, 'XeLiF')

        self.name = name 
        self.segment_area = self.hexagon_area(self.pm['segmentation_parameters']['segment_size'][0] / 2. * u.m) 
        self.total_collecting_area = self.segment_area * self.pm['segmentation_parameters']['number_segments'][0] 
        self.effective_aperture = 2. * (self.total_collecting_area / np.pi )**0.5 





