#!/usr/bin/env python
"""
Created on Fri Oct 14 20:28:51 2016
@authors: gkanarek, tumlinson
"""
import os, yaml
import math
from collections import defaultdict

from syotools.models.base import PersistentModel
from syotools.defaults import default_telescope
from syotools.spectra import utils
#from syotools.utils import pre_encode
#from syotools.utils.jsonunit import str_jsunit
import astropy.units as u #for unit conversions
import numpy as np
import scipy as sc
import synphot as syn
from hwome.core.navigator import DataModel
from syotools.models.camera import Camera
from syotools.models.ifs import IFS
from syotools.models.spectrograph import Spectrograph

class Telescope(PersistentModel):
    """
    The basic telescope class, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        name - The name of the telescope (string)
        effective_aperture - The size of the primary telescope aperture, in meters (float)
            note: there is no such thing as "aperture", there is only "effective aperture".
                For a circular/keystone primary, this is just the diameter of the circle.
                for a hex-pattern segmented primary, this is the diameter of a circle
                with the same area as the summed area of all the hex segments.
                all code should use ONLY effective_aperture
        unobscured_fraction - The fraction of the primary mirror which is not obscured (float)
        temperature - instrument temperature, in Kelvin (float)
        ota_emissivity - emissivity factor for a TMA (float)
        diff_limit_wavelength - diffraction limit wavelength, in nm (float)

        _default_model - used by PersistentModel

        cameras - the Camera objects for this telescope
    """

    def __init__(self, **kw):

        self.instruments = {}

        self.name = ''
        self.aperture = 0. * u.m
        self.temperature = 0. * u.K
        self.ota_emissivity = 0. * u.dimensionless_unscaled
        self.diff_limit_wavelength = 0. * u.nm
        self.unobscured_fraction = 1. * u.dimensionless_unscaled

        self.verbose = False
        super().__init__(default_model=default_telescope, **kw)



    # @property
    # def effective_aperture(self):
    #     unobscured, aper = self.recover('unobscured_fraction', 'aperture')
    #     return np.sqrt(unobscured) * aper

    def add_instrument(self, instrument):
        self.instruments[instrument.name] = instrument
        instrument.telescope = self

    def hexagon_area(self, side):
        return 3. * 3.**0.5 / 2. * side**2

    def set_from_sei(self, name):
        if name in ("EAC1", "EAC2", "EAC3", "EAC5"):
            tel = self.set_from_hwome(name.lower())
        else:
            print('We do not have SEI information for: ', name)
            raise NotImplementedError

    def set_from_hwome(self,name):
        self.name = name
        self.hwo_data = DataModel()
        self.hwo_data.load_hardware(f"{name}.yaml")

        self.telescope_bands = {}

        for instrument in self.hwo_data.Instrument:
            if "Coronagraph" not in instrument.name.value:
                try:
                    modenames = list(instrument.Channel.name.keys())
                except (KeyError, TypeError):
                    modenames = [f"{instrument.name.value}.HRI_A_VIS"]
                for modename in modenames:
                    if "IFU" in modename or "IFS" in modename:
                        tel_instrument = IFS(self)
                        tel_instrument.set_from_hwome(modename, "ifs")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_IFS"] = tel_instrument
                            self.telescope_bands[f"{modename}_IFS"] = tel_instrument.configuration["band"]
                    else:
                        tel_instrument = Camera(self)
                        tel_instrument.set_from_hwome(modename, "imager")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_Imager"] = tel_instrument
                            self.telescope_bands[f"{modename}_Imager"] = tel_instrument.configuration["band"]
                        tel_instrument = Spectrograph(self)
                        tel_instrument.set_from_hwome(modename, "spectrograph")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_Spectrograph"] = tel_instrument
                            self.telescope_bands[f"{modename}_Spectrograph"] = tel_instrument.configuration["band"]


        # this also sets self.effective_area
        self.effective_diameter = self.hwo_data.OTA.circumscribing_diameter.q

    def save_to_dict(self):
        output = {}
        for instrument in self.instruments:
            output[instrument] = self.instruments[instrument].save_to_dict()
        output["name"] = self.name
        output["effective_diameter"] = self.effective_diameter

        output = simplify_data(output)

        return output

    def load_from_dict(self, config):
        """
        Restore a telescope from a stored dictionary
        """
        config = complexify_data(config)

        self.name = config.pop("name")
        self.effective_diameter = config.pop("effective_diameter")

        self.instruments = {}

        for instrument in config:
            if config[instrument]["ins_type"] == "imager":
                inst = Camera(self)
            elif config[instrument]["ins_type"] == "spectrograph":
                inst = Spectrograph(self)
            elif config[instrument]["ins_type"] == "ifs":
                inst = IFS(self)
            inst.load_from_dictionary(config[instrument])
            self.instruments[instrument] = inst

    @property
    def effective_area(self):
        return self._effective_area

    @effective_area.setter
    def effective_area(self, new_area):
        # trap any values that aren't float- or float-compatible or the correct unit
        try:
            new_area/(2 * u.cm**2)
        except Exception as err:
            raise err
        if isinstance(new_area, (int, float)):
            new_area = float(new_area) << u.cm**2
        # linking them like this should ensure we always get consistent numbers
        self._effective_area = new_area
        self._effective_diameter = (np.sqrt(new_area / np.pi) * 2.).to(u.m)

    @property
    def effective_diameter(self):
        return self._effective_diameter

    @effective_diameter.setter
    def effective_diameter(self, new_diameter):
        # trap any values that aren't float- or float-compatible or the correct unit
        try:
            new_diameter/(2 * u.m)
        except Exception as err:
            raise err
        if isinstance(new_diameter, (int, float)):
            new_diameter = float(new_diameter) << u.m
        self._effective_diameter = new_diameter
        self._effective_area = (np.pi * (new_diameter/2.)**2).to(u.cm**2)

    def find_instrument_with(self, kind, wavelength=None):
        """
        Convenience function to find an instrument with specific wavelength coverage

        Parameters
        ----------
        kind : str
            "filter" or "disperser", as desired.
        wavelength : float, optional
            specific wavelength to search for, by default None

        Returns
        -------
        suitable_instruments: dict
            A dictionary of instruments, each with their list of suitable bands
        suitable_bands: dict
            A dictionary of suitable bands, each value is the instrument
        """
        suitable_instruments = defaultdict(list)
        suitable_bands = {}
        for insname in self.telescope_bands:
            for band in self.telescope_bands[insname]:
                item = self.telescope_bands[insname][band]
                if item["kind"] == kind.lower():
                    if wavelength is not None:
                        if (wavelength >= item["wave_min"]) and (wavelength <= item["wave_max"]):
                            suitable_bands[band] = insname
                            suitable_instruments[insname].append(band)
                    else:
                        suitable_bands[band] = insname
                        suitable_instruments[insname].append(band)

        return suitable_instruments, suitable_bands

    def set_from_json(self,name):
        if self.verbose:
            print('Setting Telescope to: ', name)

        if ('EAC1' in name): tel = read_json.eac1()
        if ('EAC2' in name): tel = read_json.eac2()
        if ('EAC3' in name): tel = read_json.eac3()

        self.name = tel['name']
        self.effective_aperture = tel['aperture_od'] * u.m
        self.temperature = tel['temperature_K'] * u.K
        self.diff_limited_wavelength = tel['diff_limited_wavelength'] * u.nm
        self.unobscured_fraction = tel['unobscured_fraction']

    def set_from_yaml(self, name):

        tel = read_yaml.read_hwo(name.lower())

        # the "tel" dictionary returned by read_yaml is nested, and therefore awkward
        # when summoning individual entries. And often, we do not need the individual
        # mirrors. So, we are going to break this dictionary up and carry the
        # mirrors and other pieces separately:
        self.mirrors = {}
        self.mirrors['PM'] = utils.set_coating(tel['PM']) # primary
        self.mirrors['SM'] = utils.set_coating(tel['SM']) # secondary
        self.mirrors['M3'] = utils.set_coating(tel['M3']) # tertiary
        self.mirrors['M4'] = utils.set_coating(tel['M4']) # fold mirror (?)

        # we are saving the integrated SpectralElement as an attribute
        self.telescope_efficiency = utils.mirror_efficiency(self.mirrors)

        self.pm = tel['PM']

        self.name = name
        if ('hex' in self.pm['segmentation']): # do this only if we have a hex segmented mirror
            self.segment_area = self.hexagon_area(self.pm['segmentation_parameters']['segment_size'][0] / 2. * u.m)
            self.total_collecting_area = self.segment_area * self.pm['segmentation_parameters']['number_segments'][0]
        else:
            self.total_collecting_area = np.pi * (self.pm['circumscribing_diameter'][0]/2.*u.m)**2
        self.effective_aperture = 2. * (self.total_collecting_area / np.pi )**0.5

        #WARNING!!! as of Oct 2024, the SEI database lists the diff limited wavelength
        # as a property of the camera, not the telescope. This is being set here, arbitarily, until that is fixed.
        self.diff_limited_wavelength = 0.5 * u.nm

        self.unobscured_fraction = (1. - self.pm['obscuration_ratio'][0]) * u.dimensionless_unscaled
