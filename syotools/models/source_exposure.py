#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:31:11 2017
@author: gkanarek, jt
"""
import copy
import numpy as np

import astropy.units as u
import astropy.constants as const
import scipy as sc

import synphot as syn
from synphot.models import Empirical1D, ConstFlux1D
import stsynphot as stsyn

from syotools.models.base import PersistentModel
from syotools.models.background import calc_zodi_flux
from syotools.models.profile import Profile

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
        self.wave = [np.zeros(1, dtype=float) * u.AA]
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
        if isinstance(quant, (int, float)):
            q = np.full(nb, val)
        if isinstance(quant, (u.Quantity)) and quant.isscalar:
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

    def _ensure_quantity(self, quant, unit, nb=None):
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
        quant = self._ensure_array(quant, nb=nb)
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

    def sn_box(self, band):
        """
        Function to set the percentage of flux going through an SN box of various sizes
        Used for extended sources

        Returns
        -------
        through_aperture : float
            The fraction of the total source flux going through the aperture
        aperture_pixels : float
            The number of pixels in the aperture (for correcting other properties)
        """

        wavelen = band["effective_wavelength"]
        geometry = self.source.geometry
        shape = geometry.get("geometry", "point")

        profile = Profile(self.telescope, self.instrument, geometry, wavelen)
        
        geometry_creator = {"point": profile.point_profile, "gaussian2d": profile.gaussian_profile, 
                            "sersic": profile.sersic_profile, "sersic_scale": profile.sersic_scale_profile,
                            "flat": profile.flat_profile, "power": profile.power_profile}


        x_rot, y_rot, x, y, xsamp, ysamp = profile.generate_profile()
        profile = geometry_creator[shape]()

        # now the extraction mask
        mask = self.instrument.extraction_mask(x,y, band)

        return np.sum(mask*profile), np.sum(mask)* u.pix**2


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
        1. The source 
        2. Sky background (assumed uniform across the aperture) 
        3. Thermal self-emission (assumed uniform across the aperture). 
           At the moment we only model the heat of the detector itself

        4. Dark current (assumed uniform across the aperture).
        This is the additional current flowing regardless of photons hitting the detector.
        It doesn't care about the detector QE or filter wheel.

        5. Read noise (assumed uniform across the aperture)
        The previous terms were all signal that accumulates with time. Read noise is the
        uncertainty introduced by the detector readout process itself; a fixed value per
        exposure.

        Once we've computed all of these values, we can proceed to the exposure
        time/SNR/magnitude calculations.

        At that point, the difference between imaging and spectroscopy matter.
        
        * All non-spatially-uniform components have their flux adjusted for the amount of
          the source's flux that passes through the extraction aperture (as computed by
          SourceExposure.sn_box)
        * All spatially uniform components have their flux adjusted for the number of
          pixels in the extraction box (as computed by SourceExposure.sn_box)
        """

        configuration, c_thermal, transform_flux = self.recover("instrument.configuration", "instrument._c_thermal", "instrument.transform_flux")
        pixel_scale = configuration["pixel_scale"]
        for detector in configuration["detector"]:
            dark_current = configuration["detector"]["dark_current"]
            qe = configuration["detector"]["total_qe"]
            read_noise = configuration["detector"]["read_noise"]

        if band["kind"] in ("disperser"):
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

        # set up an appropriately sized aperture
        encircled_energy, sn_box = self.sn_box(band)

        #sn_box = np.median(sn_box)

        # fsource is:
        # shaped
        # goes through the full optical path + QE
        # accumulates over time
        flux_source = source.sed * encircled_energy

        self.sky = calc_zodi_flux(source, wave, sn_box, pixel_scale)


        # fsky is:
        # uniform
        # goes through the full optical path QE
        # accumulates over time
        # Synphot doesn't like dividing a spectrum by an area unit. 
        # Rest assured, sky was supposed to be in ABMag/arcsec**2, so 
        # ABMag/arcsec**2 * pixels**2 * arcsec**2/pixel**2 is flux.
        flux_sky = self.sky * (sn_box * pixel_scale**2).value
        #print("Skyflux", flux_sky(flux_sky.waveset))

        # thermal is:
        # uniform
        # goes through the filter wheel and QE
        # accumulates over time
        thermal = c_thermal(wave, sn_box)


        #flux_source_before = sc.integrate.simpson(flux_source(flux_source.waveset), flux_source.waveset)

        # apply internal effects within telescope & instrument
        fsource = syn.observation.Observation(flux_source, band["bandpass"] * qe, binset=wave, force="taper")
        fsky = syn.observation.Observation(flux_sky, band["bandpass"] * qe, binset=wave, force="taper")
        self.thermal = syn.observation.Observation(thermal, band["bandpass"] * qe, binset=wave, force="taper")

        #flux_source_after = sc.integrate.simpson(fsource(fsource.waveset), fsource.waveset)

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
            wave = band["bandpass"].pivot()

        return wave, fsource_countrate, fsky_countrate, thermal_countrate, dark, read_noise

    def calculate(self, custom_band=None):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude
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
        Calculate for exposure times. If a custom_band has been passed in, use that.
        Otherwise, use all the bands in the channel.

        Parameters
        ----------
        custom_band : str
            Name of a band. Defaults to none.
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
        self.wave = []
        _initial_band = self.instrument.band
        _snr_temp = self._ensure_quantity(self._snr, u.dimensionless_unscaled, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            wave, result = self._update_exptime(self.source, configuration["bands"][band])
            self._exptime.append(result)
            self.wave.append(wave)
        self._snr = _snr_temp
        self.instrument.band = _initial_band

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a custom_band has been passed in, use that.
        Otherwise, use all the bands in the channel.

        Parameters
        ----------
        custom_band : str
            Name of a band. Defaults to none.
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
        self.wave = []
        _initial_band = self.instrument.band
        _exptime_temp = self._ensure_quantity(self._exptime, u.s, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self.instrument.band = band
            wave, result = self._update_snr(self.source, configuration["bands"][band])
            self._snr.append(result)
            self.wave.append(wave)
        self._exptime = _exptime_temp
        self.instrument.band = _initial_band

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Calculate for magnitudes. If a custom_band has been passed in, use that.
        Otherwise, use all the bands in the channel.

        Parameters
        ----------
        custom_band : str
            Name of a band. Defaults to none.
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
        self.wave = []
        _initial_band = self.instrument.band
        _exptime_temp = self._ensure_quantity(self._exptime, u.s, len(bands))
        _snr_temp = self._ensure_quantity(self._snr, u.dimensionless_unscaled, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            # The analytic solution is having problems right now
            # result = self._update_magnitude(self.source, configuration["bands"][band])
            wave, result = self._do_update_magnitude(self.source, configuration["bands"][band])
            self._magnitude.append(result)
            self.wave.append(wave)

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
        wave, fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        snr2 = -(_snr**2)

        a = (fsource_countrate)**2
        b = snr2 * (fsource_countrate + (fsky_countrate + thermal_countrate + dark_current)) * u.ct
        c = snr2 * read_noise**2 * _nexp
        texp = ((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).to(u.s)

        if self.verbose:
            print("Fstar:", fsource_countrate)
            print("Texp:", texp)


        exptime = texp

        return wave, exptime

    def _update_magnitude(self, source, band):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        As of 2026-09-08, this does not work as reliably as it should. It still
        produces usually-decent first guesses.
        """
        self.instrument._print_initcon(self.verbose)

        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')
        effective_area = self.recover("telescope.effective_area")
        configuration, ab_zeropoint = self.recover("instrument.configuration", "instrument.ab_zeropoint")
        qe = configuration["detector"]["total_qe"]

        # all of these are now rates, in the extraction aperture (except read_noise)
        wave, fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        # print("Fsource", fsource_countrate)


        read_noise /= u.ct**0.5
        _exptime = _exptime.to(u.s)

        snr2 = -(_snr ** 2)
        f0 = ab_zeropoint(band)
        #5509900. * (u.photon / u.s / u.cm**2) / band["bandpass"].pivot().to_value(u.nm)
        eff = 1/(band["bandpass"]*qe).efficiency()
        #eff = 1

        # bandwave = (band["bandpass"]*qe).waveset
        # bandpass = (band["bandpass"]*qe)(bandwave)
        # effband = sc.integrate.simpson(bandpass, bandwave)
        # fullband = sc.integrate.simpson(np.ones_like(bandpass), bandwave)
        # eff = effband/fullband

        a0 = (eff * _exptime)**2
        b0 = snr2 * eff * _exptime
        c0 = snr2 * ((eff * fsky_countrate + thermal_countrate + dark_current) * _exptime + (read_noise**2 * _nexp)) / u.ct
        k = (-b0 + np.sqrt(b0**2 - 4. * a0 * c0)) / (2. * a0)

        obs = syn.observation.Observation(source.sed, band["bandpass"] * qe, force="taper")

        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / obs.effective_wavelength().to(u.cm)

        #flux = k * phot_energy / effective_area

        flux = (4. * k) / (f0 * effective_area * (band["bandpass"]*qe).equivwidth().to(u.AA))
        # #flux = flux.value
        # #flux /= (band["bandpass"]*qe).efficiency()# /  3.3244442805918006
        # print("Flux/initial", fsource_countrate / (k * u.ct))
        # print("Flux/initial", fsource_countrate / flux)

        # #flux *= 45.08873273293901

        # bandwave = (band["bandpass"]*qe).waveset
        # bandpass = (band["bandpass"]*qe)(bandwave)
        # effband = sc.integrate.simpson(bandpass, bandwave)
        # fullband = sc.integrate.simpson(np.ones_like(bandpass), bandwave)
        # print((band["bandpass"]*qe).rectwidth().to(u.nm))
        # print((band["bandpass"]*qe).equivwidth().to(u.nm))
                

        # print(effband/fullband, (band["bandpass"]*qe).efficiency())
        # #sourceElement = syn.spectrum.SpectralElement(Empirical1D, points=source.sed.waveset, lookup_table=source.sed(source.sed.waveset)/np.max(source.sed(source.sed.waveset)))
        # #flux /= effband/fullband

        # #from matplotlib import pyplot as plt
        # #plt.plot((band["bandpass"]*qe).waveset, (band["bandpass"]*qe)((band["bandpass"]*qe).waveset))
        # #plt.show()

        magnitude = -2.5 * np.log10(np.array(flux)) * u.mag('AB')

        # print("eff", eff)
        # print("readnoise", read_noise**2)
        # print("dark_current", dark_current)
        # print("fsky", fsky_countrate)
        # print("Flux", flux)
        # print("Flux FNU", syn.units.convert_flux((band["bandpass"]*qe).pivot(), _magnitude, syn.units.PHOTNU))
        # print("K", k)
        # print("SNR", _snr)
        # print("A0:", a0)
        # print("B0:", b0)
        # print("C0:", c0)
        # print("F0:", f0)
        # print("Mag:", _magnitude)

        return wave, magnitude

    def _do_update_magnitude(self, source, band):
        """
        This stopgap calc-for-magnitude works differently: It sets up a range of magnitudes and modifies the source for each one.

        This is obviously super slow, as it requires computing a grid.

        Parameters
        ----------
        source : Source
            A configured source intended to be used in the calculation
        band : dict
            A bandpass dictionary
        """
        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')

        temp_magnitudes = []
        temp_snrs = []
        # make a grid of potential magnitudes covering a nice wide range
        wave, _magnitude = self._update_magnitude(source, band)
        _magnitude = _magnitude.to_value(u.ABmag)
        for temp_magnitude in np.linspace(_magnitude+4, _magnitude-2, 15):
            sp_norm = source.sed.normalize(temp_magnitude * u.ABmag, stsyn.spectrum.band(source.renorm_band))
            
            source.sed = sp_norm
            dummy, temp_snr = self._update_snr(source, band)
            temp_snrs.append(temp_snr)
            temp_magnitudes.append(temp_magnitude)
        
        maginterp = sc.interpolate.make_interp_spline(temp_snrs, temp_magnitudes, k=3)

        magnitude = maginterp(_snr) * u.ABmag

        return wave, magnitude

    def _update_snr(self, source, band):
        """
        Calculate the SNR for the given exposure time and source SED.
        """

        self.instrument._print_initcon(self.verbose)

        (_exptime, _nexp) = self.recover('_exptime', 'n_exp')

        # all of these are now rates, in the extraction aperture (except read_noise)
        wave, fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        # print("Fsource", fsource_countrate)
        # print("Fsky", fsky_countrate)
        # print("Thermal", thermal_countrate)
        # print("Dark", dark_current)
        # print("Readnoise", read_noise)

        time_per_exposure = _exptime / _nexp

        signal_counts = (fsource_countrate * _exptime).to(u.ct)
        shot_noise_in_signal = np.sqrt(signal_counts)

        sky_counts = (fsky_countrate * _exptime).to(u.ct)
        shot_noise_in_sky = np.sqrt(sky_counts)

        read_counts = (read_noise**2 * _nexp / u.ct).to(u.ct)

        dark_counts = (dark_current * _exptime).to(u.ct)

        thermal_counts = (thermal_countrate * _exptime).to(u.ct)

        tsnr = signal_counts / np.sqrt(signal_counts + sky_counts + read_counts
                                      + dark_counts + thermal_counts)
        snr = tsnr.value * u.dimensionless_unscaled

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
            
        return wave, snr

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

class SourceMultiSpecExposure(SourceExposure):
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
        Calculate for exposure times. If a custom_band has been passed in, use that. 
        Otherwise, use all of the bands in the channel.


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
        self.wave = []
        self.waves = []
        _snr_temp = self._ensure_array(self._snr, len(bands))
        # IFS and MOS instruments are valuable because they can observe multiple sources simultaneously.
        for idx,band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            _single_exptime = []
            _single_exptimemax = []
            _single_wave = []
            for source in self.sources:
                self._snr = _snr_temp[idx]
                wave, result = self._update_exptime(source, configuration["bands"][band])
                _single_exptime.append(result)
                _single_exptimemax.append(np.max(result).to_value(u.s))  # boil it down to a single scalar number to avoid numpy inhomogenous array issues
                _single_wave.append(wave)
            self._exptimes.append(_single_exptime)
            self.waves.append(_single_wave)
            # find the highest exposure time amongst the set of sources
            maxidx = np.argmax(_single_exptimemax)
            self._exptime.append(_single_exptime[maxidx])
            self.wave.append(_single_wave[maxidx])

        self._snr = _snr_temp

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a custom_band has been passed in, use that. 
        Otherwise, use all of the bands in the channel.

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
        self.wave = []
        self.waves = []
        _exptime_temp =  self._ensure_array(self._exptime, len(bands))
        # IFS and MOS instruments are valuable because they can observe multiple sources simultaneously.
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            _single_snr = []
            _single_snrmax = []
            _single_wave = []
            for source in self.sources:
                self._exptime = _exptime_temp[idx]
                wave, result = self._update_snr(source, configuration["bands"][band])
                _single_snr.append(result)
                _single_snrmax.append(np.max(result)) # boil it down to a single number to avoid numpy inhomogenous array issues
                _single_wave.append(wave)
            self._snrs.append(_single_snr)
            self.waves.append(_single_wave)
            # find the highest exposure time amongst the set of sources
            maxidx = np.argmax(_single_snrmax)
            self._snr.append(_single_snr[maxidx])
            self.wave.append(_single_wave[maxidx])

        self._exptime = _exptime_temp

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Not supported, make this an error
        """
        raise ValueError("Magnitude calculation not supported for MultiSpec Spectroscopy")

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
