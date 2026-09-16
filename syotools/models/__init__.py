#!/usr/bin/env python
"""
Created on Fri Oct 14 20:28:51 2016

@author: gkanarek, tumlinson
"""

from .telescope import Telescope
from .camera import Camera
from .spectrograph import Spectrograph, Spectropolarimeter
from .multispec import MultiSpec
from .multispec import MultiSpec as IFS
from .coronagraph import Coronagraph
#from .exposure import PhotometricExposure, SpectrographicExposure
from .source_exposure import SourcePhotometricExposure, SourceSpectrographicExposure, SourceCoronagraphicExposure, SourceMultiSpecExposure
from .source_exposure import SourceMultiSpecExposure as SourceIFSExposure
from .source import Source
