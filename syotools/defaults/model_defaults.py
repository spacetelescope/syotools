#!/usr/bin/env python
"""
Created on Sat Oct 15 10:59:16 2016

@author: gkanarek, tumlinson
"""
import os
from collections import OrderedDict
from syotools.utils import ordered_load 
import astropy.units as u 

cwd = os.getenv('LUVOIR_SIMTOOLS_DIR')

# Load data from fits files describring LUMOS and establish the default file path
current_dir = os.path.dirname(os.path.abspath(__file__))
spec_default_path = os.path.join(current_dir, '..', 'reference_data', 'LUMOS_ETC.fits')
pol_default_path = os.path.join(current_dir, '..', 'reference_data', 'POLLUX_ETC.fits')
yaml_default_path = os.path.join(current_dir, 'model_defaults.yaml')

# Load the defaults from model_defaults.yaml
#Default exposure parameters are taken from the default HDI_ETC tool values
#   --> the _sed, _snr, and _magnitude default values are placeholders
#Default spectrograph parameters are taken from the default LUMOS_ETC tool values
#   --> the description, bef, R, wrange, wave, and aeff default values are placeholders
with open(yaml_default_path, 'r') as stream:
    all_defaults = ordered_load(stream)
default_telescope = all_defaults['Telescope']
default_camera = all_defaults['Camera']
default_exposure = all_defaults['Exposure']
default_spectrograph = OrderedDict([("_lumos_default_file", spec_default_path)])
default_spectrograph.update(all_defaults['Spectrograph'])
default_spectropolarimeter = OrderedDict([("_lumos_default_file", pol_default_path)])
default_spectropolarimeter.update(all_defaults['Spectropolarimeter'])
default_coronagraph = all_defaults['Coronagraph'] #placeholder
