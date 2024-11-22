#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Th Nov 21 2024 JT 
"""
import numpy as np
import astropy.units as u

from syotools.models.base import PersistentModel
from syotools.spectra import SpectralLibrary
from syotools.spectra.utils import renorm_sed
import pysynphot as S 

class Source(PersistentModel):
    def __init__(self):
        """
        Initialize a Source object. 

        Parameters:
        - name (str): The name of the source.
        - magnitude (float): The magnitude of the source.
        - redshift (float): The redshift value of the source.
        - extinction (float): The extinction value of the source.

        By default, this object is initialized with a pysynphot 
        flat spectrum in AB mag, normalized to ABmag = 30.  

        usage: 
            > from syotools.models.source import Source
            > s = Source()

            s.sed can then be manipulated with pysynphot syntax like so: 
            > s.sed.renorm(20., 'abmag', S.ObsBandpass('johnson,v'))
        """
   
        self.name = 'fab' 
        self.magnitude = 30. 
        self.redshift = 0. 
        self.extinction = 0.  

        self.sed = S.FlatSpectrum(30, fluxunits='abmag')


    def __repr__(self):
        """
        Provide a string representation of the Source object.
        """
        return (f"Source(name={self.name!r}, magnitude={self.magnitude}, "
                f"redshift={self.redshift}, extinction={self.extinction})")




