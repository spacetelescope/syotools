#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Th Nov 21 2024 JT 
"""
import numpy as np
import astropy.units as u
import astropy.constants as const

from syotools.models.base import PersistentModel
from syotools.utils import pre_encode, pre_decode
from syotools.spectra import SpectralLibrary
from syotools.spectra.utils import renorm_sed

def nice_print(arr):
    """
    Utility to make the verbose output more readable.
    """
    
    arr = pre_decode(arr) #in case it's a JsonUnit serialization

    if isinstance(arr, u.Quantity):
        l = ['{:.2f}'.format(i) for i in arr.value]
    else:
        l = ['{:.2f}'.format(i) for i in arr]
    return ', '.join(l)

class Source(PersistentModel):
    def __init__(self):
        """
        Initialize a Source object.

        This source will be taken from the SYOTools SpectralLibrary
        which is produced by ~/spectra/library.py 

        Parameters:
        - name (str): The name of the source.
        - magnitude (float): The magnitude of the source.
        - redshift (float): The redshift value of the source.
        - extinction (float): The extinction value of the source.
        """
   
        self.sed = SpectralLibrary._available_spectra['fab']
        self.name = 'fab' 
        self.magnitude = 30. 
        self.redshift = 0. 
        self.extinction = 0.  

    def __repr__(self):
        """
        Provide a string representation of the Source object.
        """
        return (f"Source(name={self.name!r}, magnitude={self.magnitude}, "
                f"redshift={self.redshift}, extinction={self.extinction})")



