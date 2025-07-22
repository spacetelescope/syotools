#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:04:24 2017

@author: gkanarek
"""
import astropy.units as u
import pysynphot as pys
from pathlib import Path
import astropy.io.ascii as asc

#Define a new unit for spectral flux density
flambda = u.def_unit(["flambda","flam"], (u.photon / u.s / u.cm**2 / u.nm))

def renorm_sed(sed, new_magnitude, bandpass='johnson,v', waveunits='nm', fluxunits='abmag'):
    """
    Utility to renormalize an SED to a new manitude.
    """
    band = pys.ObsBandpass(bandpass)
    band.convert(sed.waveunits)
        
    new_sed = sed.renorm((new_magnitude + 2.5*u.mag('AB')).value, 'abmag', band)
    new_sed.convert(waveunits) 
    new_sed.convert(fluxunits) 
        
    return new_sed

def mag_from_sed(sed, camera):
    """
    Given a spectral energy distribution (SED) and a camera, convolve the
    SED with the camera's bandpass(es) to generate a set of AB magnitudes.
    
    NOTE AS OF 2017-11-20: not going to convolve or anything atm, just spit
    out pysynphot samples at the band pivotwaves. Convolution & integration
    over the camera bandpasses will be implemented in a future version.
    
    Parameters:
        sed    - pysynphot spectrum
        camera - Camera model
    """
    
    #Acquire camera bandpasses, making use of astropy.modeling model sets.
    pivots = camera.recover('pivotwave')
    
    sed.convert('ABMag')
    sed.convert(pivots.unit.name)
    output_mag = sed.sample(pivots.value)
    
    return output_mag * u.ABmag

# utility for when we need to load a text file  
def load_txtfile(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    path = Path(fname[0])
    for f in fname[1:]:
        path = path / f
    abspath = str(path.resolve())

    tab = asc.read(abspath, names=['wave','flux']) 
    sp = pys.ArraySpectrum(wave=tab['wave'].value, flux=tab['flux'].value, waveunits='Angstrom', fluxunits='flam')
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')
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
    sp = pys.ArraySpectrum(wave=tab['lam'].value, flux=tab['lh1=17.5'].value, 
                           waveunits='Angstrom', fluxunits='flam')
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')
    sp.__setattr__('band', band)
    return sp

# utility for when we need to load a pysynphot spectrum  
def load_pysfits(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    path = Path(fname[0])
    for f in fname[1:]:
        path = path / f
    abspath = str(path.resolve())
    sp = pys.FileSpectrum(abspath)
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')
    sp.__setattr__('band', band)
    return sp
