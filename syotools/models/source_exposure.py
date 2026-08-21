#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:31:11 2017
@author: gkanarek, jt
"""
import numpy as np
import astropy.units as u
import astropy.constants as const

import synphot as syn
from synphot.models import Empirical1D

from syotools.models.base import PersistentModel

from syotools.defaults import default_exposure
from syotools.models.source import Source

class SourceExposure(PersistentModel):
    """
    The base Source exposure class, which provides parameter storage for
    optimization, and all exposure-specific calculations. Since the
    Nov 2024 refactor, this class uses the Source object to handle
    astrophysical source information. Also, all JSON encoding has been
    stripped out.

    The SNR, exptime, and limiting magnitude can each be calculated from the
    other two. To trigger such calculations when parameters are updated, we
    will need to create property setters.

    Attributes:
        telescope    - the Telescope model instance associated with this exposure
        camera       - the Camera model instance associated with this exposure
        spectrograph - the Spectrograph model instance (if applicable) associated
                       with this exposure
        ifs          - the IFS model instance (if applicable) associated with this exposure

        exp_id       - a unique exposure ID, used for save/load purposes (string)
                        NOTE: THIS HAS NO DEFAULT, A NEW EXP_ID IS CREATED
                        WHENEVER A NEW CALCULATION IS SAVED.
        source       - Source object that this Exposure will observe
        n_exp        - the desired number of exposures (integer)
        exptime      - the desired exposure time (float array)
        snr          - the desired S/N ratio (float array)
        magnitude    - either the input source magnitude, in which case this is
                       equal to the SED interpolated to the desired wavelengths,
                       or the limiting magnitude of the exposure (float array)
        unknown      - a flag to indicate which variable should be calculated
                       ('snr', 'exptime', or 'magnitude'). this should generally
                       be set by the tool, and not be available to users. (string)
        sources      - list of source objects to be added to this exposure

        _default_model - used by PersistentModel
    """

    def __init__(self, default_model=default_exposure, **kw):

        self.source = Source() # this is the Source object, returns a flat spectrum by default.
                        # currently an Exposure can have only one Source

        self.telescope = None
        self.instrument = None

        self.exp_id = ''
        self.n_exp = 1
        self._exptime = np.ones(1, dtype=float) * u.h
        self._snr = np.zeros(1, dtype=float)
        self._magnitude = np.zeros(1, dtype=float) * u.ABmag
        self._unknown = "" # one of 'snr', 'magnitude', 'exptime'
        self._interp_flux = np.zeros(1, dtype=float) * u.dimensionless_unscaled # the source SED interpolated to the Spectrograph wavelength grid

        self.verbose = False # set this to True for debugging purposes
        self._disable = True #set this to disable recalculating (when updating several attributes at the same time)
        #super().__init__(default_model, **kw)

    def disable(self):
        self._disable = True

    def enable(self):
        self._disable = False
        self.calculate()

    #Property wrappers for the three possible unknowns, so that we can auto-
    #calculate whenever they're set, and to prevent overwriting previous
    #calculations by accident.

    @property
    def unknown(self):
        return self._unknown

    @unknown.setter
    def unknown(self, new_unknown):
        valid_unknowns = ("exptime", "snr", "magnitude")
        if new_unknown in valid_unknowns:
            self._unknown = new_unknown
            self.enable() # once this is set, enable calculation (which immediately runs one)
        else:
            raise KeyError(f"Cannot solve for {new_unknown}, unrecognized unknown.")

    def _ensure_array(self, quant, nb=None):
        """
        Ensure that the given Quantity is an array, propagating if necessary.
        """
        if self.instrument is None:
            nb = 1
        elif nb is None:
            nb = self.recover('instrument.n_bands')
        val = quant 
        if quant.isscalar:
            q = np.full(nb, val)
        elif len(quant) < nb:
            if len(quant) > 1:
                q = np.full(nb, val[0])
            else:
                q = np.full(nb, val)
        elif len(quant) > nb:
            q = quant[0:nb]
        else:
            q = val

        if not isinstance(q, u.Quantity):
            q = q << quant.unit

        return q

    def _ensure_quantity(self, quant, unit):
        """
        Ensure given quantity is an astropy unit.Quantity
        of appropriate type
        """
        if isinstance(quant, u.Quantity):
            # just see if this crashes.
            try:
                quant.to(unit)
            except:
                raise ValueError(f"Quantity {quant} unit is not convertible to {unit}.")
        else:
            quant = quant << unit
        quant = self._ensure_array(quant)
        return quant

    @property
    def exptime(self):
        return self._exptime

    @exptime.setter
    def exptime(self, new_exptime):
        if self.unknown == "exptime":
            return
        self._exptime = self._ensure_quantity(new_exptime, u.s)
        self.calculate()

    @property
    def snr(self):
        return self._snr

    @snr.setter
    def snr(self, new_snr):
        if self.unknown == "snr":
            return
        self._snr = self._ensure_quantity(new_snr, u.dimensionless_unscaled)
        self.calculate()

    @property 
    def magnitude(self):
        return self._magnitude

    @magnitude.setter 
    def magnitude(self, new_magnitude):
        if self.unknown == "magnitude":
            return
        self._magnitude = self._ensure_quantity(new_magnitude, u.ABmag)
        self.calculate()


    @property
    def interpolated_sed(self):
        """
        The exposure's (old style) SED interpolated at the camera bandpasses.
        """
        if not self.camera:
            return self.sed
        sed = self.recover('sed')
        return self.camera.interpolate_at_bands(sed)

    def process_observation(self, source, band, verbose=False):
        """
        Process the entire observation, up through the point we compute SNR/Exptime/Mag
        
        The components of flux in the observation are: 
        1. The source (assumed to be a point, but let's give it a size 
           parameter) 
        2. Sky background (assumed uniform across the aperture) 
        3. Thermal self-emission (assumed uniform across the aperture).
        At the moment we only model the heat of the detector itself

        4. Dark current (assumed uniform across the aperture).
        This is the additional current flowing regardless of photons hitting
        the detector. It doesn't care about the detector QE or filter wheel.

        5. Read noise (assumed uniform across the aperture)
        The previous terms were all signal that accumulates with time. Read
        noise is the uncertainty introduced by the detector readout process 
        itself; a fixed value per exposure.

        Once we've computed all of these values, we can proceed to the
        exposure time/SNR/magnitude calculations.

        At that point, the difference between imaging and spectroscopy matter.
        
        For imaging:
        * All non-spatially-uniform components have (size * psf size) 
        compared to (aperture size), light losses adjusted accordingly, and 
        integrated over the bandpass + QE (source, sky) or QE (thermal) to be 
        single value(s)
        * All uniform components processed for the aperture size

        For spectroscopy:
        * All non-spatially-uniform components are convolved with a 
        response function equal to the resolving power of the instrument
        and then have their (size * psf size) compared to slit size 
        (width * height, if applicable), light losses adjusted accordingly, 
        and convolved with the bandpass+QE (source, sky) or QE (thermal), 
        then convolved with a response function equal to the resolving power 
        of the instrument.
        * All uniform components processed for the height of the slit * 
        resolving power.
        """

        sky = self.recover("instrument.sky")
        configuration, c_thermal, _sn_box, transform_flux = self.recover("instrument.configuration", "instrument._c_thermal", "instrument._sn_box", "instrument.transform_flux")
        pixel_scale = configuration["pixel_scale"]
        for detector in configuration["detector"]:
            dark_current = configuration["detector"]["dark_current"]
            qe = configuration["detector"]["total_qe"]
            read_noise = configuration["detector"]["read_noise"]

        if band["kind"] in ("disperser", "ifs"):
            R = band["resolution"]
            waveunit = band["bandpass"].waveset.unit
            wavepix = np.linspace(band["bandpass"].waveset[0], band["bandpass"].waveset[-1], 1000) # using the bandpass wavelengths leads to weird fringing
            delta_lambda = wavepix/R
            pixel = np.cumsum(1.0 / delta_lambda * np.gradient(wavepix))
            pixel_integer = np.arange(int(pixel[0]), int(pixel[-1]))
            wave = np.interp(pixel_integer, pixel, wavepix) << waveunit
            # Or just use the instrument bandpass
            #wave = band["bandpass"].waveset

            dw = wave[1:] - wave[:-1]
            good = np.where(dw != 0)[0]
            wave = wave[good]
            dw = dw[good]
            #dw = np.append(dw, dw[-1])
        else:
            dw = 1
            wave = source.sed.waveset
        syn.utils.validate_wavelengths(wave)
        self.wave = wave

        # set up an appropriately sized aperture
        sn_box = _sn_box(self.wave, False)

        sn_box = np.median(sn_box)

        # fsource is:
        # shaped
        # goes through the full optical path + QE
        # accumulates over time
        flux_source = source.sed
        # scale source radius to the aperture size - we get all of the flux if it's smaller than the aperture
        if source.radius > 0:
            area = np.pi * (source.radius/pixel_scale)**2
        else:
            area = np.pi * (np.median(self.instrument.fwhm_psf(self.wave))/pixel_scale)**2
        if area > sn_box:
            flux_source = source.sed * (sn_box/area)

        # fsky is:
        # uniform
        # goes through the full optical path QE
        # accumulates over time
        # Synphot doesn't like dividing a spectrum by an area unit. 
        # Rest assured, sky was supposed to be in ABMag/arcsec**2, so 
        # ABMag/arcsec**2 * pixels**2 * arcsec**2/pixel**2 is flux.
        flux_sky = sky * (sn_box * pixel_scale**2).value
        #print("Skyflux", flux_sky(flux_sky.waveset))

        # thermal is:
        # uniform
        # goes through the filter wheel and QE
        # accumulates over time
        thermal = c_thermal(self.wave)

        # print("Source", flux_source.waveset)
        # print("Sky", flux_sky.waveset)
        # print("Thermal", thermal.waveset)
        # total_band = band["bandpass"] * qe
        # total_flux = total_band(total_band.waveset)
        # b1, b2 = total_band.waveset.min(), total_band.waveset.max()
        # a1, a2 = flux_source.waveset.min(), flux_source.waveset.max()
        # print("WAVE EDGES", a1, b1, a2, b2, a2 < b1, b2 < a1)
        # print("Disjoint", total_band.check_overlap(flux_source))
        # print("Valid band", total_band.waveset[total_flux > 0])
        # print("Band", (band["bandpass"]* qe).waveset)
        # print("QE", qe.waveset)



        # apply internal effects within telescope & instrument
        fsource = syn.observation.Observation(flux_source, band["bandpass"] * qe, binset=self.wave, force="taper")
        fsky = syn.observation.Observation(flux_sky, band["bandpass"] * qe, binset=self.wave, force="taper")
        self.thermal = syn.observation.Observation(thermal, band["bandpass"] * qe, binset=self.wave, force="taper")

        # dark is:
        # uniform
        # only within detector
        # accumulates over time
        dark = dark_current * sn_box# * qy

        # readnoise is:
        # uniform
        # only within detector
        # single event at read time
        read_noise = read_noise * np.sqrt(sn_box) # * u.electron**0.5 / u.pix**0.5

        fsource_countrate = transform_flux(fsource, wave) * dw
        fsky_countrate = transform_flux(fsky, wave) * dw
        thermal_countrate = transform_flux(self.thermal, wave) * dw
        if dw == 1:
            self.wave = band["bandpass"].pivot()

        return fsource_countrate, fsky_countrate, thermal_countrate, dark, read_noise

    def calculate(self, custom_band=None):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        """
        if self._disable:
            return False
        if self.instrument is None or self.telescope is None:
            return False
        result = {'magnitude': self.calculate_magnitude,
                'exptime': self.calculate_exptime,
                'snr': self.calculate_snr}[self.unknown](custom_band=custom_band)

        return result

    def calculate_exptime(self, custom_band=None):
        """
        Calculate for exposure times. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._exptime = []
        _initial_band = self.instrument.band
        _snr_temp = self._ensure_array(self._snr, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            result = self._update_exptime(self.source, configuration["band"][band])
            self._exptime.append(result)
        self._snr = _snr_temp
        self.instrument.band = _initial_band

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._snr = []
        _initial_band = self.instrument.band
        _exptime_temp = self._ensure_array(self._exptime, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self.instrument.band = band
            result = self._update_snr(self.source, configuration["band"][band])
            self._snr.append(result)
        self._exptime = _exptime_temp
        self.instrument.band = _initial_band

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Calculate for magnitudes. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._magnitude = []
        _initial_band = self.instrument.band
        _exptime_temp = self.exptime
        _snr_temp = self._snr
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            result = self._update_magnitude(self.source, configuration["band"][band])
            self._magnitude.append(result)

        self._exptime = _exptime_temp
        self._snr = _snr_temp
        self.instrument.band = _initial_band

        return True

    def _update_exptime(self, source, band):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        self.instrument._print_initcon(self.verbose)

        (_snr, _nexp) = self.recover('_snr', 'n_exp')

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        snr2 = -(_snr**2)

        a = (fsource_countrate)**2
        b = snr2 * (fsource_countrate + (fsky_countrate + thermal_countrate + dark_current)) * u.ct
        c = snr2 * read_noise**2 * _nexp
        texp = ((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).to(u.s)

        if self.verbose:
            print("Fstar:", fsource_countrate)
            print("Texp:", texp)


        _exptime = texp

        return _exptime

    def _update_magnitude(self, source, band):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """
        self.instrument._print_initcon(self.verbose)

        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')
        effective_area = self.recover("telescope.effective_area")
        configuration = self.recover("instrument.configuration")
        qe = configuration["detector"]["total_qe"]

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        read_noise /= u.ct**0.5
        _exptime = _exptime.to(u.s)

        snr2 = -(_snr ** 2)
        f0 = 5509900. * (u.photon / u.s / u.cm**2) / band["bandpass"].pivot()

        a0 = (_exptime)**2
        b0 = snr2 * _exptime
        c0 = snr2 * ((fsky_countrate + thermal_countrate + dark_current) * _exptime + (read_noise**2 * _nexp)) / u.ct
        k = (-b0 + np.sqrt(b0**2 - 4. * a0 * c0)) / (2. * a0)

        flux = (4. * k) / (f0 * effective_area * (band["bandpass"]*qe).equivwidth().to(u.nm))

        flux *= band["bandpass"].tlambda()

        _magnitude = -2.5 * np.log10(np.array(flux)) * u.mag('AB')

        return _magnitude

    def _update_snr(self, source, band):
        """
        Calculate the SNR for the given exposure time and source SED.
        """

        self.instrument._print_initcon(self.verbose)

        (_exptime, _nexp) = self.recover('_exptime', 'n_exp')

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        time_per_exposure = _exptime / _nexp

        signal_counts = fsource_countrate * _exptime
        shot_noise_in_signal = np.sqrt(signal_counts)

        sky_counts = fsky_countrate * _exptime
        shot_noise_in_sky = np.sqrt(sky_counts)

        read_counts = read_noise**2 * _nexp / u.ct

        dark_counts = dark_current * _exptime

        thermal_counts = thermal_countrate * _exptime

        snr = signal_counts / np.sqrt(signal_counts + sky_counts + read_counts
                                      + dark_counts + thermal_counts)
        _snr = snr.value * u.dimensionless_unscaled

        if self.verbose:
            print('# of exposures: {}'.format(_nexp))
            print('Time per exposure: {}'.format(time_per_exposure))
            print('Signal counts: {}'.format(self.nice_print(signal_counts)))
            print('Signal shot noise: {}'.format(self.nice_print(shot_noise_in_signal)))
            print('Sky counts: {}'.format(self.nice_print(sky_counts)))
            print('Sky shot noise: {}'.format(self.nice_print(shot_noise_in_sky)))
            print('Total read noise: {}'.format(self.nice_print(read_counts)))
            print('Dark current noise: {}'.format(self.nice_print(dark_counts)))
            print('Thermal counts: {}'.format(self.nice_print(thermal_counts)))
            print('SNR: {}'.format(snr))
            
        return _snr

    def add_source(self, new_source):
        self.source = new_source

class SourcePhotometricExposure(SourceExposure):
    """ A subclass of the base Exposure model, for photometric ETC calculations """
    pass



class SourceSpectrographicExposure(SourceExposure):
    """
    A subclass of the base Exposure model, for spectroscopic ETC calculations.
    """

    def calculate_magnitude(self, custom_band=None):
        """
        Not supported, make this an error
        """
        raise ValueError("Magnitude calculation not supported for Spectroscopy")

class SourceIFSExposure(SourceExposure):
    """ 
    This is currently a subclass of Spectrographic exposure that accepts multiple
    sources and produces multiple returns. 
    """
    def __init__(self, default_model=default_exposure, **kw):

        # need this before so super().__init__ has somewhere to put the default source
        self.sources = []
        self._snrs = []
        self._exptimes = []
        self._wavelength = []

        super().__init__(default_model, **kw)
        # Do this after, because by default super().__init__ loads a default source
        self.sources = []

    def add_source(self, source):
        # and now the magic: create a master wavelength array from all of the sources.
        self.sources.append(source)
        for source in self.sources:
            self._wavelength = syn.utils.merge_wavelengths(self._wavelength, syn.models.get_waveset(source.sed.model))

    @property
    def num_sources(self):
        return len(self.sources)

    @property
    def source(self):
        return self.sources[-1]

    @source.setter
    def source(self, new_source):
        self.add_source(new_source)

    @property
    def exptimes(self):
        return self._exptimes

    @exptimes.setter
    def exptimes(self, new_exptime):
        print("Did not set exptimes")

    @property
    def snrs(self):
        return self._snrs

    @snrs.setter
    def snrs(self, new_snr):
        print("Did not set snrs")

    def calculate_exptime(self, custom_band=None):
        """
        Calculate for exposure times. If a band has been passed in, do that. 
        Otherwise, do all of the bands in the channel.


        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._exptime = []
        self._exptimes = []
        _snr_temp = self._ensure_array(self._snr, len(bands))
        # The unique thing about IFS is it has multiple sources
        for source in self.sources:
            _single_exptime = []
            for idx,band in enumerate(bands):
                # because a multiple-in, multiple-out is a valid use case
                self._snr = _snr_temp[idx]
                result = self._update_exptime(source, configuration["band"][band])
                _single_exptime.append(result)
            self._exptimes.append(_single_exptime)
        # find the highest exposure time amongst the set of sources
        self._exptime = np.max(self._exptimes,axis=0)
        
        self._snr = _snr_temp

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a band has been passed in, do that. 
        Otherwise, do all of the bands in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._snr = []
        self._snrs = []
        _exptime_temp =  self._ensure_array(self._exptime, len(bands))
        # The unique thing about IFS is it has multiple sources
        for source in self.sources:
            _single_snr = []
            for idx, band in enumerate(bands):
                # because a multiple-in, multiple-out is a valid use case
                self._exptime = _exptime_temp[idx]
                result = self._update_snr(source, configuration["band"][band])
                _single_snr.append(result)
            self._snrs.append(_single_snr)
        # find the highest exposure time amongst the set of sources
        self._snr = np.max(self._snrs,axis=0)

        self._exptime = _exptime_temp

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Not supported, make this an error
        """
        raise ValueError("Magnitude calculation not supported for IFS Spectroscopy")

class SourceCoronagraphicExposure(SourceExposure):
    """
    A subclass of the base Exposure model, for coronagraphic imaging calculations.
    """

    def calculate(self):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        JT - THIS PART DOESNT WORK FOR CORON YET
        """
        if self._disable:
            return False
        if self.camera is None or self.telescope is None:
            return False
        status = {'magnitude': self._update_magnitude,
                  'exptime': self._update_exptime,
                  'snr': self._update_snr}[self.unknown]()
        return status

    #Calculation methods
    def _update_exptime(self):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        print("Doesn't exist yet, pull it from camera class")

        return False #completed successfully

    def _update_magnitude(self):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """

        print("Doesn't exist yet, pull it from camera class")

        return False #completed successfully

    def _update_snr(self):
        """
        Calculate the SNR for the given exposure time and planet properties.
        Follows Mennesson et al. 2024
        """

        self.instrument._print_initcon(self.verbose)

        print(' telescope inside the Coron exposure object ',
         self.telescope.effective_aperture)

        #serialize with JsonUnit for transportation
        self._snr = 10.

        return True #completed successfully
