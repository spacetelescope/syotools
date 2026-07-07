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

def nice_print(arr):
    """ Utility to make the verbose output more readable. """

    if isinstance(arr, u.Quantity):
        l = ['{:.2f}'.format(i) for i in arr.value]
        unit = str(arr.unit)
    else:
        l = ['{:.2f}'.format(i) for i in arr]
        unit = ''
    return ', '.join(l) + '  ' + unit

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
        self.camera = None
        self.spectrograph = None
        self.spectropolarimeter = None
        self.ifs = None

        self.exp_id = ''
        self.n_exp = 0
        self._exptime = np.zeros(1, dtype=float) * u.h
        self._snr = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self._magnitude = np.zeros(1, dtype=float) * u.ABmag
        self._unknown = '' # one of 'snr', 'magnitude', 'exptime'
        self._interp_flux = np.zeros(1, dtype=float) * u.dimensionless_unscaled # the source SED interpolated to the Spectrograph wavelength grid

        self.verbose = True # set this to True for debugging purposes
        self._disable = False #set this to disable recalculating (when updating several attributes at the same time)
        super().__init__(default_model, **kw)

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
        self._unknown = new_unknown
        self.calculate()

    def _ensure_array(self, quant):
        """
        Ensure that the given Quantity is an array, propagating if necessary.
        """
        q = quant
        #val = q[1]['value']
        val = q # not sure about this. - it should be stripping out the JSON but leaving the intent intact.
        if not isinstance(val, list):
            if self.camera is None:
                nb = 1
            else:
                nb = self.recover('camera.n_bands')
            q[1]['value'] = np.full(nb, val).tolist()

        return q

    @property
    def exptime(self):
        #print(" retrieve exptime")
        return self._exptime

    @exptime.setter
    def exptime(self, new_exptime):
        if self.unknown == "exptime":
            return
        self._exptime = self._ensure_array(new_exptime)
        self.calculate()

    @property
    def snr(self):
        return self._snr

    @snr.setter
    def snr(self, new_snr):
        if self.unknown == "snr":
            return
        self._snr = self._ensure_array(new_snr)
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

    def interpolated_source(self, source):
        """
        The exposure's new Source SED interpolated at the camera bandpasses.

        telescope efficiency reduces counts at detector (HWOE-183)
        """
        thru, qe = self.recover("camera.throughput", "camera.total_qe")
        output_mags = [] # <--- create blank list of mags
        for band in thru:
            # multiply the sed by the bandpass
            bandpass = band["bandpass"]
            sed = syn.observation.Observation(source.sed, bandpass)
            # extract the magnitude in AB Magnitudes
            this_mag = sed.effstim(u.ABmag)
            output_mags.append(this_mag.value)
            if self.verbose:
                print('getting mags from interpolated _source: ', bandpass.avgwave())
        return np.array(output_mags)

    def process_observation(self, source, band, verbose=True):
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
        (widwth * height, if applicable), light losses adjusted accordingly, 
        and convolved with the bandpass+QE (source, sky) or QE (thermal), 
        then convolved with a response function equal to the resolving power 
        of the instrument.
        * All uniform components processed for the height of the slit * 
        resolving power.
        """
        wave = source.sed.waveset
        # set up an appropriately sized aperture
        Npix = self.instrument._sn_box(wave, False)
        thru, qe, c_thermal, dark_current, readnoise = self.recover("instrument.throughput", "instrument.total_qe", "instrument.c_thermal", "instrument.dark_current", "instrument.readnoise")

        # fsource is:
        # shaped
        # goes through the full optical path + QE
        # accumulates over time
        flux_source = syn.units.convert_flux(source.sed.waveset, source.sed(source.sed.waveset), syn.units.PHOTLAM)
        if source.radius > 0:
            radius = source.radius
        else:
            radius = instrument.fwhm_psf(wave)
        if radius > np.sqrt(Npix)/2.:
            fsource = fsource * (np.sqrt(Npix)/2)/radius

        # fsky is:
        # uniform
        # goes through the full optical path QE
        # accumulates over time
        flux_sky = syn.units.convert_flux(sky.sed.waveset, sky.sed(sly.sed.waveset), syn.units.PHOTLAM)
        flux_sky *= Npix
        
        # thermal is:
        # uniform
        # goes through the filter wheel and QE
        # accumulates over time
        thermal = c_thermal(wave)
        Omega = (pixel_size**2 * box * u.pix).to(u.sr)
        thermal *= Omega

        # apply internal effects within telescope & instrument
        fsource = syn.observation.Observation(flux_source, band["bandpass"] * qe)
        fsky = syn.observation.Observation(flux_sky, band["bandpass"] * qe)
        thermal = syn.observation.Observation(thermal, band["bandpass"] * qe)

        # dark is:
        # uniform
        # only within detector
        # accumulates over time
        dark = dark_current

        # readnoise is:
        # uniform
        # only within detector
        # single event at read time
        readnoise = readnoise

        return fsource, fsky, thermal, dark, readnoise


    @property 
    def magnitude(self, source=None):
        if self.unknown == "magnitude":
            return self._magnitude
        #If magnitude is not unknown, it should be interpolated from the SED #at the
        camera bandpasses. if self.verbose:
            print('magnitude fcn line 191', self.interpolated_source(source))
        return self.interpolated_source(source)

    @magnitude.setter 
    def magnitude(self, new_magnitude):
        if self.unknown == "magnitude":
            return
        self._magnitude = self._ensure_array(new_magnitude) if self.verbose:
            print('magnitude fcn line 200', new_magnitude)

        self.calculate()

    def calculate(self):
        """
        This method is implemented at the subclass level, not here
        (so in PhotometricExposure, SpectroscopicExposure, etc.)
        """
        raise NotImplementedError

    def add_source(self, new_source):
        self.source = new_source

class SourcePhotometricExposure(SourceExposure):
    """ A subclass of the base Exposure model, for photometric ETC calculations """

    def calculate(self):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        """
        if self._disable:
            return False
        if self.camera is None or self.telescope is None:
            return False
        status = {'magnitude': self._update_magnitude,
                  'exptime': self._update_exptime,
                  'snr': self._update_snr}[self.unknown](self.source)
        return status

    def _fsource(self, source):
        """
        Calculate the stellar flux as per Eq 2 in the SNR equation paper.
        """
        mag = self.interpolated_source(source)
        (f0, c_ap, D, dlam) = self.recover('camera.ab_zeropoint',
                                           'camera.ap_corr',
                                           'telescope.effective_aperture',
                                           'camera.derived_bandpass')

        m = 10.**(-0.4*(mag)) # magnitude to flux
        D = D.to(u.cm)

        fsource = f0 * c_ap[0] * np.pi / 4. * D**2 * (dlam * u.nm) * m

        return fsource

    def _fsky(self, verbose=True):
        """
        Calculate the sky flux as per Eq 6 in the SNR equation paper.
        """

        (f0, D, dlam, Phi, fwhm, Sigma, throughput, qe, pivotwave) = self.recover('camera.ab_zeropoint',
                'telescope.effective_aperture', 'camera.derived_bandpass', 'camera.pixel_size', 
                'camera.fwhm_psf', 'camera.sky_sigma', 'camera.throughput', 'camera.total_qe', 'camera.pivotwave')

        D = D.to(u.cm)
        m = 10.**(-0.4 * np.array(Sigma)) / u.arcsec**2
        Npix = self.camera._sn_box(False)

        if verbose:
            print('Sky brightness: {}'.format(Sigma))

        fsky = f0 * np.pi / 4. * D**2 * (dlam*u.nm) * m * (Phi**2 * Npix) * u.pix
        # telescope efficiency reduces counts at detector (HWOE-183)

        for bidx, band in enumerate(throughput):
            bandpass = band["bandpass"]
            fsky[bidx] *= bandpass(pivotwave[bidx])

        return fsky

    def _update_exptime(self, source):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        self.camera._print_initcon(self.verbose)

        (_snr, _nexp) = self.recover('_snr', 'n_exp')
        (throughput, effective_diameter) = self.recover("instrument.throughput", "telescope.effective_aperture")
        effective_aperture = (effective_diameter/2)**2 * np.pi
        for band in throughput:
            fsource, fsky, thermal, dark, readnoise = process_observation(source, band)
            fsource = fsource.countrate(effective_aperture)
            fsky = fsky.countrate(effective_aperture)

        print("SNR in function", _snr)
        snr2 = -(_snr**2)
        fstar = self._fsource(source)
        fsky = self._fsky(verbose=self.verbose)
        Npix = self.camera._sn_box(self.verbose)
        thermal = self.camera.c_thermal(verbose=self.verbose)

        dark_rate = _dark_current * u.electron/u.count
        rn = _detector_rn.value * u.electron**0.5/u.pixel**0.5

        QE = _total_qe(pivotwave) * u.electron/u.photon
        a = (QE * fstar)**2
        b = snr2 * (QE * (fstar + fsky) + thermal + dark_rate * Npix)
        c = snr2 * rn**2 * Npix * _nexp
        texp = ((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).to(u.s)

        if self.verbose:
            print("Fstar:", fstar)
            print("Texp:", texp)

        self._exptime = texp

        return True

    def _update_magnitude(self, source):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """
        self.camera._print_initcon(self.verbose)

        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')
        (f0, c_ap, D, dlam, pivotwave) = self.recover('camera.ab_zeropoint',
                                           'camera.ap_corr',
                                           'telescope.effective_aperture',
                                           'camera.derived_bandpass',
                                           "camera.pivotwave")
        (_total_qe, _detector_rn, _dark_current, pivotwave) = self.recover('camera.total_qe',
                                                            'camera.readnoise', 
                                                            'camera.dark_current', 
                                                            "camera.pivotwave")

        exptime = (_exptime[0] * u.Unit(_exptime[1])).to(u.s)

        D = D.to(u.cm)
        fsky = self._fsky(verbose=self.verbose)

        Npix = self.camera._sn_box(self.verbose)
        c_t = self.camera.c_thermal(verbose=self.verbose)

        QE = _total_qe(pivotwave) * u.electron/u.photon
        rn = _detector_rn.value * u.electron**0.5/u.pixel**0.5
        dark_rate = _dark_current * u.electron/u.count

        snr2 = -(_snr ** 2)

        a0 = (QE * exptime)**2
        b0 = snr2 * QE * exptime
        c0 = snr2 * ((QE * fsky + c_t + Npix * dark_rate) * exptime + (rn**2 * Npix * _nexp))
        k = (-b0 + np.sqrt(b0**2 - 4. * a0 * c0)) / (2. * a0)
        # telescope efficiency reduces the flux at the detector, so it must be divided out
        # (increase flux) to find the actual flux required to get that SNR in that time.
        # (HWOE-183)
        flux = (4. * k) / (f0 * c_ap[0] * np.pi * D**2 * (dlam*u.nm))# / int_eff(pivotwave[0] * u.Unit(pivotwave[1]))

        self._magnitude = -2.5 * np.log10(np.array(flux)) * u.mag('AB')

        return True

    def _update_snr(self, source):
        """
        Calculate the SNR for the given exposure time and source SED.
        """

        self.camera._print_initcon(self.verbose)

        (_exptime, _nexp, n_bands) = self.recover('_exptime', 'n_exp',
                                                  'camera.n_bands')

        (_total_qe, _detector_rn, _dark_current, pivotwave) = self.recover('camera.total_qe',
                'camera.readnoise', 'camera.dark_current', "camera.pivotwave")

        number_of_exposures = np.full(n_bands, _nexp)
        desired_exp_time = (np.full(n_bands, _exptime[0]) * u.Unit(_exptime[1])).to(u.second)
        time_per_exposure = desired_exp_time / number_of_exposures

        QE = _total_qe(pivotwave) * u.electron/u.photon

        signal_counts = QE * self._fsource(source) * desired_exp_time
        shot_noise_in_signal = np.sqrt(signal_counts)

        # telescope efficiency reduces counts at detector (HWOE-183)
        sky_counts = QE * self._fsky(verbose=self.verbose) * desired_exp_time
        shot_noise_in_sky = np.sqrt(sky_counts)

        sn_box = self.camera._sn_box(self.verbose) #<-- units should be "pix"

        rn = _detector_rn.value * u.electron**0.5/u.pixel**0.5
        read_counts = rn**2 * sn_box * number_of_exposures

        dark_rate = _dark_current * u.electron/u.count
        dark_counts = sn_box * dark_rate * desired_exp_time

        thermal_counts = desired_exp_time * self.camera.c_thermal(verbose=self.verbose)

        snr = signal_counts / np.sqrt(signal_counts + sky_counts + read_counts
                                      + dark_counts + thermal_counts)
        self._snr = snr

        if self.verbose:
            print('# of exposures: {}'.format(_nexp))
            print('Time per exposure: {}'.format(time_per_exposure[0]))
            print('Signal counts: {}'.format(nice_print(signal_counts)))
            print('Signal shot noise: {}'.format(nice_print(shot_noise_in_signal)))
            print('Sky counts: {}'.format(nice_print(sky_counts)))
            print('Sky shot noise: {}'.format(nice_print(shot_noise_in_sky)))
            print('Total read noise: {}'.format(nice_print(read_counts)))
            print('Dark current noise: {}'.format(nice_print(dark_counts)))
            print('Thermal counts: {}'.format(nice_print(thermal_counts)))
            print('SNR: {}'.format(snr))
            print('Max SNR: {} in {} band'.format(snr.max(), self.camera.bandnames[snr.argmax()]))

        return True

class SourceSpectrographicExposure(SourceExposure):
    """
    A subclass of the base Exposure model, for spectroscopic ETC calculations.
    """

    def calculate(self):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        """
        if self._disable:
            return False
        if self.spectrograph is None or self.telescope is None:
            return False

        if self.unknown == "snr":
            self._update_snr(self.source)
        if self.unknown == "exptime":
            self._update_exptime(self.source)

    def _update_snr(self, source):
        """
        Calculate the SNR based on the current SED and spectrograph parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('effective_aperture'))
            msg2 = " with {} in mode {}".format(self.spectrograph.name, self.spectrograph.mode)
            print(msg1 + msg2)

        _exptime = self.recover('exptime')
        _wave, aeff, bef, aper, R, wrange = self.recover('spectrograph.wave',
                                                         'spectrograph.aeff',
                                                         'spectrograph.bef',
                                                         'telescope.effective_aperture',
                                                         'spectrograph.R',
                                                         'spectrograph.wrange')

        exptime = ( self._exptime[0][0] * u.Unit(self._exptime[1])).to(u.s)

        wave = _wave.to(u.AA)

        sourceobs = syn.observation.Observation(source.sed, aeff)
        

        sflux = syn.units.convert_flux(sourceobs.binset, sourceobs.binflux, (u.erg / u.s / u.cm**2 / u.AA))

        delta_lambda = self.recover('spectrograph.delta_lambda').to(u.AA / u.pix)
        pixel = np.cumsum(1.0 / delta_lambda * np.gradient(wave))
        pixel_integer = np.arange(int(pixel.value[0]), int(pixel.value[-1]))
        #wavepix = np.interp(pixel_integer, pixel.value, wave)
        wavepix = wave

        # with open("wavefile.csv", "w") as outfile:
        #     for wave in wavepix:
        #         outfile.write(f"{wave.value}\n")


        #iflux = np.interp(wave, swave, sflux, left=0., right=0.)
        # telescope efficiency reduces counts at detector (HWOE-183)
        #self._interp_flux = iflux
        #phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.pix

        source_counts = sourceobs.countrate(np.pi * (aper/2)**2, wavelengths=wavepix) * exptime
        #print('_update_snr phot_energy: ', phot_energy)

        # This is the aperture efficiency (throughput) scaled from 15m radius to the requested aperture size
        #scaled_aeff = aeff * (aper / (15 * u.m))**2
        #source_counts = iflux / phot_energy * scaled_aeff * exptime * delta_lambda
        #print('_update_snr source_counts: ', source_counts)

        # telescope efficiency reduces counts at detector (HWOE-183)
        bgobs = syn.observation.Observation(bef, aeff)
        bg_counts = bgobs.countrate(np.pi * (aper/2)**2, wavelengths=wavepix) * exptime
        #bg_counts = bef / phot_energy * scaled_aeff * exptime

        snr = source_counts / np.sqrt(source_counts + bg_counts)

        if self.verbose:
            print("SNR: {}".format(snr))

        self._snr = snr

        return True

    def _update_exptime(self, source):
        """
        Calculate the exptime based on the current SED and spectrograph parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('aperture'))
            msg2 = " with {} in mode {}".format(self.spectrograph.name, self.spectrograph.mode)
            print(msg1 + msg2)

        _snr, _exptime = self.recover('_snr', '_exptime')
        _wave, aeff, bef, aper, R, wrange= self.recover('spectrograph.wave',
                                                         'spectrograph.aeff',
                                                         'spectrograph.bef',
                                                         'telescope.effective_aperture',
                                                         'spectrograph.R',
                                                         'spectrograph.wrange',
                                                         )

        if self.verbose:
            print("The requested SNR is {}\n".format(_snr))

        wave = _wave.to(u.AA)
        wavepix = wave

        scaled_aper = np.pi * (aper/2)**2

        sourceobs = syn.observation.Observation(source.sed, aeff)
        source_counts = sourceobs.countrate(scaled_aper, wavelengths=wavepix)
        bgobs = syn.observation.Observation(bef, aeff)
        bg_counts = bgobs.countrate(scaled_aper, wavelengths=wavepix)

        #sflux = syn.units.convert_flux(swave, source.sed(swave), u.erg / u.s / u.cm**2 / u.AA)

        delta_lambda = self.recover('spectrograph.delta_lambda').to(u.AA / u.pix)

        #iflux = np.interp(wave, swave, sflux, left=0., right=0.)
        # telescope efficiency reduces flux at detector (HWOE-183)
        #iflux *= telescope_efficiency
        #bef *= telescope_efficiency

        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.ct

        scaled_aeff = aeff * (aper / (15 * u.m))**2

        if (self.verbose):
            print('sflux  = ', source.sed(source.sed.waveset), '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
            print('wave = ', wave, '\n') #<---- 20,600 element array of wavelengths tied to Spectrograph object (not Exposure)
            print('delta_lambda = ', delta_lambda, '\n') #<--- this has the correct units, "Angstrom/pix"
            print('iflux = ', sourceobs.binflux, '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
                                 #<--- because the units are carried through the interpolation
            print('bef = ', bef(bef.waveset))  #<--- this has the correct units, "erg / (pix s cm2)"
            print('photE = ', phot_energy, '\n') #<--- this has the correct units, "erg / ct"
            print('aeff = ', aeff) #<--- this has the correct units, "cm2"
            print('aper = ', aper)#<--- this has the correct units, "m"
            print('scaled_aeff = ', scaled_aeff, '\n') #<--- this has the correct units, "cm2"
            print('SNR^2 :', (_snr)**2)

        t_exp = (_snr)**2 * (source_counts + bg_counts) / (source_counts**2)

        if self.verbose:
            print("Exptime: {}".format(t_exp))

        self._exptime = t_exp

        return True

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

    def calculate(self):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        """
        if self._disable:
            return False
        if self.ifs is None or self.telescope is None:
            return False

        if self.unknown == "snr":
            self._update_snrs(self.source)
        if self.unknown == "exptime":
            self._update_exptimes(self.source)

    def add_source(self, source):
        # and now the magic: create a master wavelength array from all of the sources.
        self.sources.append(source)
        for source in self.sources:
            self._wavelength = syn.utils.merge_wavelengths(self._wavelength, syn.models.get_waveset(source.sed.model))

    def _update_snr(self, source):
        """
        Calculate the SNR based on the current SED and IFS parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('aperture'))
            msg2 = " with {} in mode {}".format(self.ifs.name, self.ifs.mode)
            print(msg1 + msg2)

        _exptime = self.recover('exptime')
        _wave, aeff, bef, aper, R, wrange = self.recover('ifs.wave',
                                                         'ifs.aeff',
                                                         'ifs.bef',
                                                         'telescope.effective_aperture',
                                                         'ifs.R',
                                                         'ifs.wrange')

        exptime = ( self._exptime[0][0] * u.Unit(self._exptime[1])).to(u.s)

        wave = _wave.to(u.AA)
        wavepix = wave

        sourceobs = syn.observation.Observation(source.sed, aeff)

        sflux = syn.units.convert_flux(sourceobs.binset, sourceobs.binflux, (u.erg / u.s / u.cm**2 / u.AA))

        delta_lambda = self.recover('ifs.delta_lambda').to(u.AA / u.pix)
        pixel = np.cumsum(1.0 / delta_lambda * np.gradient(wave))
        pixel_integer = np.arange(int(pixel.value[0]), int(pixel.value[-1]))

        #iflux = np.interp(wave, swave, sflux, left=0., right=0.)
        # telescope efficiency reduces counts at detector (HWOE-183)
        #iflux *= internal_efficiency
        #self._interp_flux = iflux
        #phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.pix

        source_counts = sourceobs.countrate(np.pi * (aper/2)**2, wavelengths=wavepix) * exptime
        #print('_update_snr phot_energy: ', phot_energy)

        # This is the aperture efficiency (throughput) scaled from 15m radius to the requested aperture size
        #scaled_aeff = aeff * (aper / (15 * u.m))**2
        #source_counts = iflux / phot_energy * scaled_aeff * exptime * delta_lambda
        #print('_update_snr source_counts: ', source_counts)

        # telescope efficiency reduces counts at detector (HWOE-183)
        bgobs = syn.observation.Observation(bef, aeff)
        bg_counts = bgobs.countrate(np.pi * (aper/2)**2, wavelengths=wavepix) * exptime
        #bg_counts = bef / phot_energy * scaled_aeff * exptime

        snr = source_counts / np.sqrt(source_counts + bg_counts)

        if self.verbose:
            print("SNR: {}".format(snr))

        self._snr = snr

        return True

    def _update_exptime(self, source):
        """
        Calculate the exptime based on the current SED and IFS parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('aperture'))
            msg2 = " with {} in mode {}".format(self.ifs.name, self.ifs.mode)
            print(msg1 + msg2)

        _snr, _exptime = self.recover('_snr', '_exptime')
        #_snr = _snr[0][0]
        _wave, aeff, bef, aper, R, wrange = self.recover('ifs.wave',
                                                         'ifs.aeff',
                                                         'ifs.bef',
                                                         'telescope.effective_aperture',
                                                         'ifs.R',
                                                         'ifs.wrange')

        if self.verbose:
            print("The requested SNR is {}\n".format(_snr))

        wave = _wave.to(u.AA)
        wavepix = wave

        scaled_aper = np.pi * (aper/2)**2

        sourceobs = syn.observation.Observation(source.sed, aeff)
        source_counts = sourceobs.countrate(scaled_aper, wavelengths=wavepix)
        bgobs = syn.observation.Observation(bef, aeff)
        bg_counts = bgobs.countrate(scaled_aper, wavelengths=wavepix)

        #sflux = syn.units.convert_flux(swave, source.sed(swave), u.erg / u.s / u.cm**2 / u.AA)

        delta_lambda = self.recover('ifs.delta_lambda').to(u.AA / u.pix)

        #iflux = np.interp(wave, swave, sflux, left=0., right=0.)
        # telescope efficiency reduces flux at detector (HWOE-183)
        #iflux *= telescope_efficiency
        #bef *= telescope_efficiency

        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.ct

        scaled_aeff = aeff * (aper / (15 * u.m))**2

        if (self.verbose):
            print('sflux  = ', source.sed(source.sed.waveset), '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
            print('wave = ', wave, '\n') #<---- 20,600 element array of wavelengths tied to Spectrograph object (not Exposure)
            print('delta_lambda = ', delta_lambda, '\n') #<--- this has the correct units, "Angstrom/pix"
            print('iflux = ', sourceobs.binflux, '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
                                 #<--- because the units are carried through the interpolation
            print('bef = ', bef(bef.waveset))  #<--- this has the correct units, "erg / (pix s cm2)"
            print('photE = ', phot_energy, '\n') #<--- this has the correct units, "erg / ct"
            print('aeff = ', aeff(aeff.waveset)) #<--- this has the correct units, "cm2"
            print('aper = ', aper)#<--- this has the correct units, "m"
            print('scaled_aeff = ', scaled_aeff(scaled_aeff.waveset), '\n') #<--- this has the correct units, "cm2"
            print('SNR^2 :', (_snr)**2)

        t_exp = (_snr)**2 * (source_counts + bg_counts) / (source_counts**2)

        if self.verbose:
            print("Exptime: {}".format(t_exp))

        self._exptime = t_exp

        return True


    @property
    def num_sources(self):
        return len(self.sources)

    @property
    def source(self):
        return self.sources[-1]

    @source.setter
    def source(self, new_source):
        self.add_source(new_source)

    def _update_exptimes(self, source):
        self._exptimes = []
        if self.sources == []:
            self._exptimes = [None]
            self._exptime = None
        else:
            # loop through by setting all the sources.
            for source in self.sources:
                self._update_exptime(source)
                self._exptimes.append(self._exptime)
            
            # and set the regular exposure time to the maximum
            # the output is a spectrum, so we want to find the highest entire
            # spectrum, not any specific value.
            self._exptime = sorted(self._exptimes, key=(lambda a: np.nanmean(a)))[-1]

    def _update_snrs(self, source):
        self._snrs = []
        if self.sources == []:
            self._snrs = [None]
            self._snr = None
        else:
            # loop through by setting all the sources.
            for source in self.sources:
                self._update_snr(source)
                self._snrs.append(self._snr)
            
            # and set the regular snr to the minimum
            # the output is a spectrum, so we want to find the lowest entire
            # spectrum, not any specific value.
            self._snr = sorted(self._snrs, key=(lambda a: np.nanmean(a)))[0]

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

        self.camera._print_initcon(self.verbose)

        print(' telescope inside the Coron exposure object ',
         self.telescope.effective_aperture)

        #serialize with JsonUnit for transportation
        self._snr = 10.

        return True #completed successfully
