#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Oct 22 2024 
@author: JT
Based on exposure.py 
"""

import numpy as np
import astropy.units as u
import astropy.constants as const

from syotools.models.base import PersistentModel
from syotools.defaults import default_exposure
from syotools.utils import pre_encode, pre_decode
from syotools.spectra import SpectralLibrary
from syotools.spectra.utils import renorm_sed

    return ', '.join(l)

class Scene(PersistentModel):
    """
    The base scene class, which provides property, parameter, and 
    data storage for 2D astronomial scenes. 
    
    This class encompasses both imaging and spectral exposures -- we will
    assume that all calculations are performed with access to the full SED 
    (which is simply interpolated to the correct wavebands for imaging).
    
    The SNR, exptime, and limiting magnitude can each be calculated from the
    other two. To trigger such calculations when parameters are updated, we
    will need to create property setters.
    
    Attributes:
        nx_pix       - number of pixels in the first dimension 
        ny_pix       - number of pixels in the second dimension 
        band         - the redshift of the SED (float)
        flux         - flux per bin/pixel in physical units 
        
        _default_model - used by PersistentModel
    """
    
    _default_model = default_exposure
    
    _nx_pix = 10 
    _ny_pix = 10 

    _magnitude = pre_encode(np.zeros([_nx_pix, _ny_pix]), dtype=float) * u.ABmag)
    _flambda   = pre_encode(np.zeros([_nx_pix, _ny_pix]), dtype=float) * u.erg / u.cm**2 / u.s / u.AA) 
 
    verbose = False #set this for debugging purposes only
    _disable = False #set this to disable recalculating (when updating several attributes at the same time)

# below here, code is just stolen from Exposure class 
    
    def _ensure_array(self, quant):
        """
        Ensure that the given Quantity is an array, propagating if necessary.
        """
        q = pre_encode(quant)
        if len(q) < 2:
            import pdb; pdb.set_trace()
        val = q[1]['value']
        if not isinstance(val, list):
            if self.camera is None:
                nb = 1
            else:
                nb = self.recover('camera.n_bands')
            q[1]['value'] = np.full(nb, val).tolist()
        
        return q
    
    def sed(self):
        """
        Return a spectrum, redshifted if necessary. We don't just store the 
        redshifted spectrum because pysynphot doesn't save the original, it 
        just returns a new copy of the spectrum with redshifted wavelengths.
        """
        sed = pre_decode(self._sed)
        z = self.recover('redshift')
        return pre_encode(sed.redshift(z))
    
    @sed.setter
    def sed(self, new_sed):
        self._sed = pre_encode(new_sed)
        self.calculate()
        
    
    @sed_id.setter
    def sed_id(self, new_sed_id):
        if new_sed_id == self._sed_id:
            return
        self._sed_id = new_sed_id
        self._sed = pre_encode(SpectralLibrary.get(new_sed_id, SpectralLibrary.fab))
        self.calculate()
        
    def renorm_sed(self, new_mag, bandpass='johnson,v'):
        sed = self.recover('_sed')
        self._sed = renorm_sed(sed, pre_decode(new_mag), bandpass=bandpass)
        self.calculate()
    
    @property
    def interpolated_sed(self):
        """
        The exposure's SED interpolated at the camera bandpasses.
        """
        if not self.camera:
            return self.sed
        sed = self.recover('sed')
        return pre_encode(self.camera.interpolate_at_bands(sed))
    
    @property
    def magnitude(self):
        if self.unknown == "magnitude":
            return self._magnitude
        #If magnitude is not unknown, it should be interpolated from the SED
        #at the camera bandpasses. 
        return self.interpolated_sed
    
    @magnitude.setter
    def magnitude(self, new_magnitude):
        if self.unknown == "magnitude":
            return
        self._magnitude = self._ensure_array(new_magnitude)
        self.calculate()
