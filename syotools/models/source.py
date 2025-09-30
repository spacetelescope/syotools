#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Th Nov 21 2024 JT 
"""
import numpy as np
import astropy.units as u

from syotools.models.base import PersistentModel
from syotools.spectra.spec_defaults import syn_spectra_library 
import synphot as syn
import stsynphot as stsyn

class Source(PersistentModel):
    def __init__(self):
        """
        Initialize a Source object. 

        Parameters:
        - name (str): The name of the source.
        - magnitude (float): The magnitude of the source.
        - redshift (float): The redshift value of the source.
        - extinction (float): The extinction value of the source.

        By default, this object is initialized with a synphot 
        flat spectrum in AB mag, normalized to ABmag = 30.  

        usage: 
            > from syotools.models.source import Source
            > s = Source()
            > s.set_sed('Flat (AB)', 25., 0., 0.0, 'galex,fuv')   
            or 
            > s.set_sed('QSO', 25., 0.0, 0.0, 'galex,fuv')   

            s.sed can also be manipulated with synphot syntax like so: 
            > s.sed.renorm(20., 'abmag', S.ObsBandpass('johnson,v'))
        """
        self.name = 'Flat (AB)'
        self.magnitude = 30. 
        self.redshift = 0. 
        self.extinction = 0.  
        self.renorm_band = 'johnson,v'

        #set default here
        self.sed = None # Will be set in set_sed, do this so sed is in __init__.
        self.set_sed(self.name, self.magnitude, self.redshift, self.extinction)

        # yes this is weird, required because the superclass expects
        # attributes to be present and initialized.
        super().__init__()


    def set_sed(self, source_name, magnitude, redshift, extinction, bandpass=None):   
        self.name = source_name  
        self.sed = syn_spectra_library[source_name]
        self.magnitude = magnitude
        self.redshift = redshift
        self.extinction = extinction
        # if the bandpass is none/unspecified, load the library default
        if bandpass is None:
            self.renorm_band = syn_spectra_library[source_name].band
        else:
            self.renorm_band = bandpass

        #print("SET SED:", bandpass, syn_spectra_library[source_name].band, self.renorm_band, stsyn.band(self.renorm_band).waveset)
        #print("SED_INFO:", self.name, self.sed.waveset, self.renorm_band, self.redshift, self.extinction)

        new_sed = syn_spectra_library[source_name]

        # now apply the other quantities via synphot 
        new_sed.z = self.redshift
        sp_ext = new_sed * syn.reddening.ReddeningLaw.from_extinction_model('mwavg').extinction_curve(self.extinction)

        #print("Actual norm:", sp_ext.waveset)

        sp_norm = sp_ext.normalize(self.magnitude * u.ABmag, stsyn.spectrum.band(self.renorm_band))
        

        self.sed = sp_norm

    def list_templates(self): 
        print(syn_spectra_library.keys()) 

    def __repr__(self):
        """
        Provide a string representation of the Source object.
        """
        return (f"Source(name={self.name!r}, magnitude={self.magnitude}, "
                f"redshift={self.redshift}, extinction={self.extinction}, "
                f"renorm band={self.renorm_band}) ")
