#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:04:24 2017

@author: gkanarek
"""
import astropy.units as u
import synphot as syn
import stsynphot as stsyn
from pathlib import Path
import astropy.io.ascii as asc

#Define a new unit for spectral flux density
flambda = u.def_unit(["flambda","flam"], (u.photon / u.s / u.cm**2 / u.nm))

def renorm_sed(sed, new_magnitude, bandpass='johnson,v', waveunits='nm', fluxunits='abmag'):
    """
    Utility to renormalize an SED to a new manitude.
    """
    band = stsyn.band(bandpass)
        
    new_sed = sed.normalize(new_magnitude * u.Unit(fluxunits),  band=band)
        
    return new_sed

def mag_from_sed(sed, camera):
    """
    Given a spectral energy distribution (SED) and a camera, convolve the
    SED with the camera's bandpass(es) to generate a set of AB magnitudes.
    
    NOTE AS OF 2017-11-20: not going to convolve or anything atm, just spit
    out synphot samples at the band pivotwaves. Convolution & integration
    over the camera bandpasses will be implemented in a future version.
    
    Parameters:
        sed    - synphot spectrum
        camera - Camera model
    """
    
    #Acquire camera bandpasses, making use of astropy.modeling model sets.
    pivots = camera.recover('pivotwave')
    
    return sed(pivots).to_value(u.ABmag)

# utility for when we need to load a text file  
def load_txtfile(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    path = Path(fname[0])
    for f in fname[1:]:
        path = path / f
    abspath = str(path.resolve())

    sp = syn.spectrum.SourceSpectrum.from_file(abspath)
    sp = sp.normalize(30.0 * u.ABmag, stsyn.band(band))
    sp.__setattr__('band', band)
    return sp 

# utility for when we need to load an fesc data file  
def load_fesc(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')

    path = Path(fname[0])
    for f in fname[1:]:
        path = path / f
    abspath = str(path.resolve())
    tab = asc.read(abspath)
    sp = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=tab['lam'], lookup_table=tab['lh1=17.5'])

    sp = sp.normalize(30. * u.ABmag, stsyn.band(band))
    sp.__setattr__('band', band)
    return sp

# utility for when we need to load a synphot spectrum  
def load_synfits(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    path = Path(fname[0])
    for f in fname[1:]:
        path = path / f
    abspath = str(path.resolve())
    sp = syn.spectrum.SourceSpectrum.from_file(abspath)
    sp = sp.normalize(30.*u.ABmag, stsyn.band(band))
    sp.__setattr__('band', band)
    return sp
