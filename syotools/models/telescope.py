#!/usr/bin/env python
"""
Created on Fri Oct 14 20:28:51 2016
@authors: gkanarek, tumlinson
"""
import os, yaml
import copy
import math
from collections import defaultdict
from importlib import metadata
from importlib import metadata
import subprocess

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
from syotools.models.multispec import MultiSpec
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
        self.name = name.lower()
        self.hwo_data = DataModel()
        self.hwo_data.load_hardware(f"{self.name}.yaml")

        self.telescope_bands = {}

        for instrument in self.hwo_data.Instrument:
            if "Coronagraph" not in instrument.name.value:
                try:
                    modenames = list(instrument.Channel.name.keys())
                except (KeyError, TypeError):
                    modenames = [f"{instrument.name.value}.HRI_A_VIS"]
                for modename in modenames:
                    if "IFU" in modename.upper() or "IFS" in modename.upper():
                        tel_instrument = MultiSpec(self)
                        tel_instrument.set_from_hwome(modename, "ifs")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_IFS"] = tel_instrument
                            self.telescope_bands[f"{modename}_IFS"] = tel_instrument.bands
                    elif "MOS" in modename.upper():
                        tel_instrument = MultiSpec(self)
                        tel_instrument.set_from_hwome(modename, "mos")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_MOS"] = tel_instrument
                            self.telescope_bands[f"{modename}_MOS"] = tel_instrument.bands
                    else:
                        tel_instrument = Camera(self)
                        tel_instrument.set_from_hwome(modename, "imager")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_Imager"] = tel_instrument
                            self.telescope_bands[f"{modename}_Imager"] = tel_instrument.bands
                        tel_instrument = Spectrograph(self)
                        tel_instrument.set_from_hwome(modename, "spectrograph")
                        if tel_instrument.configuration["channel_filters"] != []:
                            self.instruments[f"{modename}_Spectrograph"] = tel_instrument
                            self.telescope_bands[f"{modename}_Spectrograph"] = tel_instrument.bands


        #print(self.hwo_data.OTA.circumscribing_diameter.q)
        #self.effective_diameter = self.hwo_data.OTA.circumscribing_diameter.q
        # This value is backed by a function that computes whether the primary mirror is
        # made of hexagons, keystones, and whether it's on-axis (with a cutout) or not.
        # The effective diameter within is based off this value assuming a perfect circle.
        self.effective_area = self.hwo_data.OTA.inscribed_aperture_area.q

    def save_to_dict(self):
        output = {}
        for instrument in self.instruments:
            output[instrument] = self.instruments[instrument].save_to_dict()
        output["name"] = self.name
        output["effective_diameter"] = self.effective_diameter

        # tag the software version the dict was created with, too
        output["syotools_version"] = metadata.version('syotools')
        output["hwome_version"] = metadata.version('hwome-core')
        output["data_version"] = subprocess.run(["git", "-C", os.environ["HWOME_DATA_PATH"], "rev-parse", "HEAD"])

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
        self._effective_area = new_area.to(u.cm**2)
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

    def find_instrument_with(self, instrument=None, kind=None, wavelength=None, resolution=None):
        """
        Convenience function to find a band (and its instrument) that meets specific
        criteria.

        Parameters
        ----------
        instrument: str, optional
            Name string found in an instrument
        kind : str, optional
            "filter" or "disperser", as desired.
        wavelength : float or list, optional
            specific wavelength to search for, by default None
        resolution : float or list, optional

        Returns
        -------
        suitable_instruments: dict
            A dictionary of instruments, each with their list of suitable bands
        suitable_bands: dict
            A dictionary of suitable bands, each value is the instrument
        """

        # options = search_configuration(channel_type='mos',
        # wavelength_range_nm=[wave_minmin/10, wave_maxmax/10],
        # resolution_range = [8000, 20000],
        # center_nm=None)

        # fbyctr = {}
        # for chan_name, cdict in options.items():
        #     for filt_name, fdict in cdict.items():
        #         print(f"Found {chan_name}.{filt_name}, center={fdict['center']}, width={fdict['width']}, R={fdict['spectral_resolution']}")
        #         fbyctr[fdict['center'].to('nm').value] = fdict

        suitable_instruments = defaultdict(list)
        suitable_bands = {}

        # set up some lists for progressive filtering
        filter_list = []
        temp_filter_list = []

        # the initial sift - literally everything
        for insname in self.telescope_bands:
            for band in self.telescope_bands[insname]:
                filter_list.append([insname, band, self.telescope_bands[insname][band]])
        # filter 1: the kind of band
        if kind is not None:
            for entry in filter_list:
                insname = entry[0]
                band = entry[1]
                item = entry[2]
                if item["kind"] == kind.lower():
                    temp_filter_list.append((insname, band, item))
            filter_list = copy.deepcopy(temp_filter_list)
            temp_filter_list = []

        # filter 2: the wavelength of the band
        if wavelength is not None:
            for entry in filter_list:
                insname = entry[0]
                band = entry[1]
                item = entry[2]        
                if isinstance(wavelength, (int, float)):
                    if (wavelength * u.AA >= item["wave_min"]) and (wavelength * u.AA <= item["wave_max"]):
                        temp_filter_list.append((insname, band, item))
                elif isinstance(wavelength, (tuple, list)):
                    if (wavelength[0] * u.AA >= item["wave_min"]) and (wavelength[1] * u.AA <= item["wave_max"]):
                        temp_filter_list.append((insname, band, item))
                elif isinstance(wavelength, dict):
                    if (wavelength["wave_min"] * u.AA >= item["wave_min"]) and (wavelength["wave_max"] * u.AA <= item["wave_max"]):
                        temp_filter_list.append((insname, band, item))
            filter_list = copy.deepcopy(temp_filter_list)
            temp_filter_list = []


        # filter 3: the resolution of the band
        if resolution is not None:
            for entry in filter_list:
                insname = entry[0]
                band = entry[1]
                item = entry[2]
                if "resolution" in item:
                    if isinstance(resolution, (int, float)):
                        if (resolution <= item["resolution"]):
                            temp_filter_list.append((insname, band, item))
                    elif isinstance(resolution, (tuple, list)):
                        if (resolution[0] <= item["resolution"]) and (resolution[1] >= item["resolution"]):
                            temp_filter_list.append((insname, band, item))
                    elif isinstance(resolution, dict):
                        if (resolution["min"] <= item["resolution"]) and (resolution["max"] >= item["resolution"]):
                            temp_filter_list.append((insname, band, item))
            filter_list = temp_filter_list
            temp_filter_list = []

        # filter 4: the instrument
        if instrument is not None:
            for entry in filter_list:
                insname = entry[0]
                band = entry[1]
                item = entry[2]
                if instrument in insname:
                    temp_filter_list.append((insname, band, item))
            filter_list = temp_filter_list
            temp_filter_list = []

        for item in filter_list:
            insname = item[0]
            band = item[1]
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
