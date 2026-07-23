#!/usr/bin/env python
"""
Created on Fri Oct 14 21:31:18 2016
@author: gkanarek, tumlinson
"""
import numpy as np
import astropy.constants as const
import astropy.units as u
import synphot as syn
import stsynphot as stsyn
from synphot.models import Empirical1D

from .instrument import Instrument
from syotools.models.source_exposure import SourcePhotometricExposure
from syotools.defaults import default_camera
from syotools.spectra.utils import mag_from_sed, mirror_efficiency, set_coating
from hwo_sci_eng.utils import read_yaml
from hwome.core.navigator import DataModel

def nice_print(arr):
    """
    Utility to make the verbose output more readable.
    """

    if isinstance(arr, u.Quantity):
        l = ['{:.2f}'.format(i) for i in arr.value]
    else:
        l = ['{:.2f}'.format(i) for i in arr]
    return ', '.join(l)

class Camera(Instrument):
    """
    The basic camera class, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with the camera
        exposures    - the list of Exposures taken with this camera

        name         - name of the camera (string)
        n_bands      - number of wavelength bands (int)
        n_channels   - number of channels (int)
        pivotwave    - central wavelengths for bands, in nanometers (float array)
        bandnames    - names of bands (string list)
        channels     - grouping of bands into channels [UV, Optical, IR],
                       and indicating the reference band for pixel size (list of tuples)
        fiducials    - fiducial wavelength of the band, for reference (float array)
        total_qe     - total quantum efficiency in each band (float array)
        ap_corr      - magnitude correction factor for aperture size (float array)
        bandpass_r   - resolution in each bandpass (float array)
        dark_current - dark current values in each band (float array)
        detector_rn  - read noise for the detector in each band (float array)
        sky_sigma    - sky background emission (float array)

        _default_model - used by PersistentModel

    The following are attributes I haven't included, and the justification:
        R_effective - this doesn't seem to be used anywhere
    """

    def __init__(self, telescope, **kw):

        self.telescope = telescope
        self._band = None
        self.exposures = []
        self.name = ''
        #self.pivotwave = np.zeros(1, dtype=float) * u.nm
        self.bandnames = ['']
        self.channels = [([],0)]
        self.fiducials = np.zeros(1, dtype=float) * u.nm
        self.total_qe = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.ap_corr = np.ones(1, dtype=float) * u.dimensionless_unscaled
        self.bandpass_r = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.dark_current = np.zeros(1, dtype=float) * (u.electron / u.s / u.pixel)
        #self.detector_rn = np.zeros(1, dtype=float) * (u.electron / u.pixel)**0.5
        self.sky = syn.spectrum.SourceSpectrum(Empirical1D, points=[0.1,10000, 20000] << u.AA, lookup_table=[24,24,24] << u.ABmag) # Hardcode a 24th magnitude background. Should actually be ABMag/arcsec^2.
        self.sky = self.sky.normalize(24 * u.ABmag, stsyn.spectrum.band("johnson,v"))
        #super().__init__(default_camera, **kw)

    @property
    def n_bands(self):
        return len(self.configuration["channel_filters"])

    @property
    def n_channels(self):
        # this has always referred to the filters
        return len(self.configuration["channel_filters"])

    @property
    def band(self):
        return self._band

    @band.setter
    def band(self, new_band):
        if new_band in self.configuration["channel_filters"]:
            self._band = new_band
        else:
            raise KeyError(f"Cannot set band {new_band}, valid options for this channel are {self.configuration["channel_filters"]}")

    @property
    def pivotwave(self):
        pivot = []
        for bpass in self.configuration["element"]:
            band = self.configuration["element"][bpass]
            pivotval = band["bandpass"].pivot()
            pivotunit = pivotval.unit
            pivot.append(pivotval.value)

        print("Pivot", pivot)

        return np.asarray(pivot) << pivotunit


    @property
    def derived_bandpass(self):
        """
        Calculate the bandpasses.
        """
        pivotwave = self.recover('pivotwave')
        width = []
        for bpass in self.configuration["element"]:
            band = self.configuration["element"][bpass]
            widthval = band["bandpass"].equivwidth()
            widthunit = widthval.unit
            width.append(widthval.value)
        width = width << widthunit

        return np.array(pivotwave/width)

    @property
    def ab_zeropoint(self):
        """
        AB-magnitude zero points as per Marc Postman's equation.
        """
        pivotwave = self.recover('pivotwave')
        pivot = pivotwave.to(u.nm)
        abzp = 5509900. * (u.photon / u.s / u.cm**2) / pivot

        print("AB Zero", abzp[0])

        return abzp# << abunit

    def _print_initcon(self, verbose):
        if verbose: #These are our initial conditions
            print('Telescope diffraction limit: {}'.format(self.telescope.diff_limit_wavelength))
            print('Telescope effective_diameter: {}'.format(self.telescope.effective_diameter))
            print('Instrument temperature: {}'.format(self.configuration["detector"]["thermal"]))
            print('Pivot waves: {}'.format(nice_print(self.pivotwave)))
            print('Pixel sizes: {}'.format(self.configuration["pixel_scale"]))
            print('AB mag zero points: {}'.format(nice_print(self.ab_zeropoint)))
            #print('Quantum efficiency: {}'.format(nice_print(self.configuration["detector"]["total_qe"](self.pivotwave))))
            print('Aperture correction: {}'.format(self.ap_corr))
            #print('Bandpass resolution: {}'.format(nice_print(self.bandpass_r[0] * u.Unit(self.bandpass_r[1]))))
            print('Derived_bandpass: {}'.format(nice_print(self.derived_bandpass)))
            print('Detector read noise: {}'.format(self.configuration["detector"]["read_noise"]))
            print('Dark rate: {}'.format(self.configuration["detector"]["dark_current"]))

    def _sn_box(self, wave, verbose):
        """
        Calculate the number of pixels in the SNR computation box.
        """

        Phi = self.configuration["pixel_scale"]
        sn_box = np.round(3. * self.fwhm_psf(wave) / Phi)

        if verbose:
            print('PSF width: {}'.format(nice_print(self.fwhm_psf(wave))))
            print('SN box width: {}'.format(nice_print(sn_box)))

        return sn_box**2


    def interpolate_at_bands(self, sed):
        """
        Interpolate an SED to obtain magnitudes for the camera's wavebands.
        """
        return mag_from_sed(sed, self)

    def interpolate_source_at_bands(self, source):
        """
        Interpolate an SED to obtain magnitudes for the camera's wavebands.
        """
        return mag_from_source(self, source)

    def create_exposure(self, source=None):
        new_exposure = SourcePhotometricExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure

    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.instrument = self
        exposure.telescope = self.telescope
        exposure.calculate()

    def transform_flux(self, spectrum, wave):
        effective_area = self.recover("telescope.effective_area")
        return spectrum.countrate(effective_area)

    def set_from_hwome(self, channelname):
        self.configuration = {}
        try:
            instrument, channel = channelname.split(".")
            instrument_data = getattr(self.telescope.hwo_data, instrument)
        except ValueError:
            raise ValueError("Need a name + channel (e.g. 'HRI-S.NIR')")
        except KeyError:
            raise KeyError(f"Unrecognized Instrument {instrument}.\n Legal values are {self.telescope.hwo_data.Instrument.name}")

        instrument_data = self.telescope.hwo_data[instrument]

        self.name = channel
        instrument_data.Channel
        try:
            channel_data = getattr(instrument_data, channel)
        except KeyError:
            raise KeyError(f"Unrecognized Channel {channel}.\n Legal values are {instrument_data.Channel.name.keys()}")
        channel_data = instrument_data[channel]


        self.configuration["element"] = {}
        self.configuration["channel_filters"] = []
        for channel_filter in channel_data.Filter.name.keys():
            filter_name = channel_filter.split(".")[-1]
            self.configuration["channel_filters"].append(filter_name)
            # this commands HWOME to walk down the entire optical path of the telescope down to the filter(grating) and collect all of the optics.
            thru = self.telescope.hwo_data.OpticalPath.select(instrument=instrument, channel=channel, filter = filter_name).throughput(include_detector=False)
            # then we multiply all of them together
            total_throughput = np.prod(thru.q, axis=0)
            # and store for later retrieval
            band = syn.spectrum.SpectralElement(Empirical1D, points=thru.w, lookup_table=total_throughput)
            wavemin = band.avgwave() - band.rectwidth()/2
            wavemax = band.avgwave() + band.rectwidth()/2
            self.configuration["element"][filter_name] = {"name": filter_name, "bandpass": band, "effective_wavelength": band.avgwave(), "wave_min": wavemin, "wave_max": wavemax, "optics": len(thru.value.keys())}

        self.configuration["diffraction_limit"] = channel_data.diffraction_limited.q
        self.configuration["pixel_scale"] = channel_data.plate_scale.q
        self.configuration["fov_x"] = channel_data.fov_x.q
        self.configuration["fov_y"] = channel_data.fov_y.q

        self.configuration["detector"] = {}
        for detector in channel_data.Detector:
            self.configuration["detector"]["name"] = detector.name
            print("READNOISE", detector.read_noise.q)
            self.configuration["detector"]["read_noise"] = detector.read_noise.q
            self.configuration["detector"]["thermal"] = detector.temperature.q
            self.configuration["detector"]["dark_current"] = detector.dark_current.q / u.pix
            w = detector.qe.w
            t = detector.qe.q
            self.configuration["detector"]["total_qe"] = syn.spectrum.SpectralElement(Empirical1D, points=w, lookup_table=t)

    def set_to_dict(self, config):
        self.configuration = config

    def get_from_dict(self):
        return self.configuration

    def set_from_sei(self, name): 

        if ('HRI' in name): hri = read_yaml.hri()
        
        # the "hri" dictionary returned by read_yaml is nested, and therefore awkward  
        # when summoning individual entries. And often, we do not need the individual 
        # components. So, we are going to break this dictionary up and carry the 
        # pieces separately:  

        self.uvis = hri['UVIS']  

        self.nir = hri['NIR'] 

        self.uvis_mirrors = {}
        for mirror in range(1, self.uvis["n_refl_optics"][0] + 1):
            self.uvis_mirrors[f"mirror{mirror}"] = set_coating(self.uvis)

        self.instrument_efficiency_uvis = mirror_efficiency(self.uvis_mirrors)

        self.nir_mirrors = {}
        for mirror in range(1, self.nir["n_refl_optics"][0] + 1):
            self.nir_mirrors[f"mirror{mirror}"] = set_coating(self.nir)

        self.instrument_efficiency_nir = mirror_efficiency(self.nir_mirrors)

