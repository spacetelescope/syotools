#!/usr/bin/env python
"""
Created on Tue Oct 18 14:10:45 2016
@author: gkanarek, jt
"""
import os

import numpy as np
import synphot as syn
import stsynphot as stsyn
import astropy.units as u
import astropy.io.ascii as asc

from syotools.utils import pre_encode
from syotools.spectra.utils import load_txtfile, load_fesc, load_synfits

# path names to find reference files 
pysyn_path = os.environ['PYSYN_CDBS']
data_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','reference_data'))

# this is where the spectral library is compiled
# these dictionary entries specify the template name, filename, band, etc. 
# then we will create two dictionaries - one containing the JSON-encoded 
# spectral templates ("default_spectra") and one containing the 
# synphot objects ("syn_spectra_library"), eventually the former 
# will be deprecated when the JSON encoding is removed. 

specs = {'Classical T Tauri': {'desc': 'Classical T-Tauri Star', 
                  'file': [data_path, 'CTTS_etc_d140pc_101116.txt'], 
                  'band': 'galex,fuv'},
         'ctts2': {'desc': 'Classical T-Tauri Star', 
                  'file': [data_path, 'ttauri.txt'],
                  'band': 'galex,fuv'},
         'M1 Dwarf': {'desc': 'M1 Dwarf', 
                    'file': [data_path, 'dM1_etc_d5pc_101116.txt'], 
                    'band': 'galex,fuv'},
         'mdwarf2': {'desc': 'M Dwarf', 
                    'file': [data_path, 'dM_d1pc_pollux.txt'],
                    'band': 'galex,fuv'}, 
         'G Dwarf': {'desc': 'G Dwarf',
                    'file': [data_path, 'dG_d5pc_pollux.txt'],
                    'band': 'galex,fuv'},
         'O5V Star': {'desc': 'O5V Star', 
                 'file': [pysyn_path, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_1.fits']},
         'G2V Star': {'desc': 'G2V Star', 
                 'file': [pysyn_path, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_26.fits']},
         'B5V Star': {'desc': 'B5V Star', 
                 'file': [pysyn_path, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_6.fits']},
         'M2V Star': {'desc': 'M2V Star', 
                 'file': [pysyn_path, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_40.fits']},
         'G191B2B (WD)': {'desc': 'G191B2B (WD)', 
                     'file': [pysyn_path, 'calspec', 'g191b2b_mod_010.fits'], 
                     'band': 'galex,fuv'},
         'GD71 (WD)': {'desc': 'GD71 (WD)',
                  'file': [pysyn_path, 'calspec', 'gd71_fos_003.fits'],
                  'band': 'galex,fuv'},
         'GD153 (WD)': {'desc': 'GD153 (WD)',
                   'file': [pysyn_path, 'calspec', 'gd153_fos_003.fits'],
                   'band': 'galex,fuv'},
         '10 Myr Starburst': {'desc': '10 Myr Starburst', 
                 'file': [data_path, '10Myr_Starburst_nodust.dat'],
                 'band': 'galex,fuv'},
         'QSO': {'desc': 'QSO', 'file': [pysyn_path, 'grid', 'agn', 'qso_template.fits'], 
                 'band': 'galex,fuv'},
         'Seyfert 1': {'desc': 'Seyfert 1', 'file': [data_path, 'Seyfert_1_template.txt']},
         'Seyfert 2': {'desc': 'Seyfert 2', 'file': [data_path,'Seyfert_2_template.txt']},
         'Liner': {'desc': 'Liner', 'file': [data_path, 'LINER_template.txt']},
         'Orion Nebula': {'desc': 'Orion Nebula', 
                   'file': [pysyn_path, 'grid', 'galactic', 'orion_template.fits']},
         'Starburst, No Dust': {'desc': 'Starburst, No Dust', 
                    'file': [pysyn_path, 'grid', 'kc96', 'starb1_template.fits']},
         'Starburst, E(B-V) = 0.6': {'desc': 'Starburst, E(B-V) = 0.6', 
                  'file': [pysyn_path, 'grid', 'kc96', 'starb6_template.fits']},
         'Elliptical Galaxy': {'desc': 'Elliptical Galaxy', 'file': [pysyn_path, 'grid', 'kc96', 
                                'elliptical_template.fits']},
         'Sbc Galaxy': {'desc': 'Sbc Galaxy', 'file': [pysyn_path, 'grid', 'kc96', 'sc_template.fits']},
         'NGC 1068': {'desc': 'NGC 1068',
                     'file': [pysyn_path, 'grid', 'agn', 'ngc1068_template.fits']},
         'Galaxy with f_esc, HI=1, HeI=1': {'desc': 'Galaxy with f_esc, HI=1, HeI=1',
                     'file': [data_path, 'fesc', 'fe_lyccontrot1.000hi1.000hei.txt'],
                     'band': 'galex,fuv'},
         'Galaxy with f_esc, HI=0.001, HeI=1': {'desc': 'Galaxy with f_esc, HI=0.001, HeI=1',
                     'file': [data_path, 'fesc', 'fe_lyccontrot0.001hi1.000hei.txt'],
                     'band': 'galex,fuv'}
        }

#default_spectra is the object that other routines will import to
#use the library 
default_spectra = {'specs':{}, 'descs':{}}
syn_spectra_library = {} # this will gradually become the replacement for default_spectra, with the whole thing in synphot form 

# now iterate over the dictionary above to load the spectra
# from text, fits, or fesc format depending on the type of file 
for (specid, spec) in specs.items():
    default_spectra['descs'][specid] = spec['desc']
    fesc = 'fesc' in spec['file']
    fits = spec['file'][-1].endswith('fits')
    if fits:
        default_spectra['specs'][specid] = load_synfits(spec)
    elif fesc:
        default_spectra['specs'][specid] = load_fesc(spec) 
    else:
        default_spectra['specs'][specid] = load_txtfile(spec) 

# this loop calls the same routines but gets back synphot objects instead of pre_encoded 
#so these will populate syn_spectra_library
for (specid, spec) in specs.items():
    default_spectra['descs'][specid] = spec['desc']
    fesc = 'fesc' in spec['file']
    fits = spec['file'][-1].endswith('fits')
    if fits:
        syn_spectra_library[spec['desc']] = load_synfits(spec) 
    elif fesc:
        syn_spectra_library[spec['desc']] = load_fesc(spec) 
    else:
        syn_spectra_library[spec['desc']] = load_txtfile(spec)         

wave = np.linspace(100,30000,300) << u.Angstrom

# now add a few other special cases from synphot (syn)
# and store them in the pre_encoded JSON format \
flatsp = syn.spectrum.SourceSpectrum(syn.models.ConstFlux1D, amplitude=30.0 * u.ABmag)
flatsp = flatsp.normalize(30.0 * u.ABmag, band=stsyn.band("johnson,v"))
# This library is built on the idea of storing flux and wavelength arrays, so we need to make this an empirical spectrum.
flatsp = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=flatsp(wave))
default_spectra['specs']['fab'] = pre_encode(flatsp)
default_spectra['descs']['fab'] = 'Flat (AB)'
flatsp.__setattr__('band', 'johnson,v')
syn_spectra_library['Flat (AB)'] = flatsp

flamsp = syn.spectrum.SourceSpectrum(syn.models.ConstFlux1D, amplitude=30.0 * syn.units.FLAM)
flamsp = flamsp.normalize(30.0 * u.ABmag, band=stsyn.band("johnson,v"))
# This library is built on the idea of storing flux and wavelength arrays, so we need to make this an empirical spectrum.
flamsp = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=flamsp(wave))
default_spectra['specs']['flam'] = pre_encode(flamsp)
default_spectra['descs']['flam'] = 'Flat in F_lambda'
flamsp.__setattr__('band', 'johnson,v')
syn_spectra_library['Flat in F_lambda'] = flamsp

bb = syn.spectrum.SourceSpectrum(syn.models.BlackBody1D, temperature=5000)
wave = np.linspace(100,3000,300) << u.Angstrom
bb = bb.normalize(30.0 * u.ABmag, band=stsyn.band('galex,fuv'))
# This library is built on the idea of storing flux and wavelength arrays, so we need to make this an empirical spectrum.
bb = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=bb(wave))
default_spectra['specs']['Blackbody5000'] = pre_encode(bb)
default_spectra['descs']['Blackbody5000'] = 'Blackbody (5000K)'
bb.__setattr__('band', 'galex,fuv')
syn_spectra_library['Blackbody (5000K)'] = bb

bb = syn.spectrum.SourceSpectrum(syn.models.BlackBody1D, temperature=100000)
bb = bb.normalize(30.0 * u.ABmag, band=stsyn.band('galex,fuv'))
# This library is built on the idea of storing flux and wavelength arrays, so we need to make this an empirical spectrum.
bb = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=bb(wave))
default_spectra['specs']['Blackbody100000'] = pre_encode(bb)
default_spectra['descs']['Blackbody100000'] = 'Blackbody (100,000K)'
bb.__setattr__('band', 'galex,fuv')
syn_spectra_library['Blackbody (100,000K)'] = bb 
