#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Th Nov 21 2024 JT 
"""
import numpy as np
import astropy.units as u
from astropy import coordinates as coord

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
        self._ra = 135.0 * u.deg
        self._dec = 20.0 * u.deg

        #set default here
        self.sed = None # Will be set in set_sed, do this so sed is in __init__.
        self.set_sed(self.name, self.magnitude, self.redshift, self.extinction)
        self.radius = 0 # point source

        # yes this is weird, required because the superclass expects
        # attributes to be present and initialized.
        super().__init__()

    @property
    def ra(self):
        return self._ra

    @ra.setter
    def ra(self, new_ra):
        if isinstance(new_ra, str):
            self._ra = coord.Angle(new_ra)
        elif isinstance(new_ra, (int, float)):
            self._ra = coord.Angle(new_ra * u.deg)
        elif isinstance(new_ra, u.Quantity):
            self._ra = coord.Angle(new_ra)
        else:
            raise ValueError(f"Unrecognized RA angle {new_ra}")

    @property
    def dec(self):
        return self._dec

    @dec.setter
    def dec(self, new_dec):
        if isinstance(new_dec, str):
            self._dec = coord.Angle(new_dec)
        elif isinstance(new_dec, (int, float)):
            self._dec = coord.Angle(new_dec * u.deg)
        elif isinstance(new_dec, u.Quantity):
            self._dec = coord.Angle(new_dec)
        else:
            raise ValueError(f"Unrecognized DEC angle {new_dec}")

    @property
    def coords(self):
        return coord.SkyCoord(ra=self.ra, dec=self.dec, frame="icrs")

    @coords.setter
    def coords(self):
        pass

    @property
    def coordinates(self):
        return self.coords

    @coordinates.setter
    def coordinates(self):
        pass

    def set_sed(self, source_name, magnitude, redshift, extinction, bandpass=None, radius=0, ra=135.0, dec=20.0, geometry={}, library=syn_spectra_library):
        self.name = source_name  
        self.sed = library[source_name]
        self.magnitude = magnitude
        self.redshift = redshift
        self.extinction = extinction
        # if the bandpass is none/unspecified, load the library default
        if bandpass is None:
            self.renorm_band = library[source_name].band
        else:
            self.renorm_band = bandpass

        # Set a radius for extended sources. 0 = unresolved point source.
        self.radius = radius

        self.ra = ra
        self.dec = dec

        # one of "point", "flat", "gaussian2d", "sersic"
        self.geometry = {"shape": "sersic", "major": 0.2 * u.arcsec, "minor": 0.1 * u.arcsec,
                        "norm_method": "surf_center", "sersic_index": 2, "surf_area_units": "arcsec^2"}
        # shape is one of "point" or "gaussian2d" or "flat"
        # norm_method is one of "surf_center", "surf_scale", "integ_infinity"
        # surf_area_units is one of "arcsec^2" or "sr"

        #print("SET SED:", bandpass, library[source_name].band, self.renorm_band, stsyn.band(self.renorm_band).waveset)
        #print("SED_INFO:", self.name, self.sed.waveset, self.renorm_band, self.redshift, self.extinction)

        new_sed = library[source_name]

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
