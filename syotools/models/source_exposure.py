#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:31:11 2017
@author: gkanarek, jt
"""
import numpy as np
import astropy.units as u
import astropy.constants as const

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

        self.exp_id = ''
        self.n_exp = 0
        self._exptime = np.zeros(1, dtype=float) * u.h
        self._snr = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self._snr_goal = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self._magnitude = np.zeros(1, dtype=float) * u.ABmag
        self._unknown = '' # one of 'snr', 'magnitude', 'exptime'
        self._interp_flux = np.zeros(1, dtype=float) * u.dimensionless_unscaled # the source SED interpolated to the Spectrograph wavelength grid

        self.verbose = False # set this to True for debugging purposes
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
        if len(q) < 2:
            import pdb; pdb.set_trace()
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

    @property
    def interpolated_source(self):
        """
        The exposure's new Source SED interpolated at the camera bandpasses.
        """
        self.source.sed.convert(self.camera.pivotwave[1]) # <---temporarily convert the sed into the units of the pivotwaves
        output_mags = [] # <--- create blank list of mags
        for magwave in self.camera.pivotwave[0]:
            this_mag = self.source.sed.sample(magwave)
            #amazingly, the sample method on the pysynphot sed does not check wavelengh limits!
            if (magwave > np.max(self.source.sed.wave)): this_mag = 99.
            if (magwave < np.min(self.source.sed.wave)): this_mag = 99.
            output_mags.append(this_mag)
            if self.verbose:
                print('getting mags from interpolated _source: ', magwave)
        return np.array(output_mags)

    @property
    def magnitude(self):
        if self.unknown == "magnitude":
            return self._magnitude
        #If magnitude is not unknown, it should be interpolated from the SED
        #at the camera bandpasses.
        if self.verbose:
            print('magnitude fcn line 174', self.interpolated_source)
        return self.interpolated_source

    @magnitude.setter
    def magnitude(self, new_magnitude):
        if self.unknown == "magnitude":
            return
        self._magnitude = self._ensure_array(new_magnitude)
        if self.verbose:
            print('magnitude fcn line 182', new_magnitude)

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
                  'snr': self._update_snr}[self.unknown]()
        return status

    @property
    def _fsource(self):
        """
        Calculate the stellar flux as per Eq 2 in the SNR equation paper.
        """
        mag = self.recover('magnitude')
        (f0, c_ap, D, dlam) = self.recover('camera.ab_zeropoint',
                                           'camera.ap_corr',
                                           'telescope.effective_aperture',
                                           'camera.derived_bandpass')

        m = 10.**(-0.4*(mag))
        D = D.to(u.cm)

        fsource = f0 * c_ap[0] * np.pi / 4. * D**2 * (dlam * u.nm) * m

        return fsource

    def _update_exptime(self):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        self.camera._print_initcon(self.verbose)

        (_snr, _nexp) = self.recover('snr', 'n_exp')
        (_total_qe, _detector_rn, _dark_current) = self.recover('camera.total_qe',
                'camera.detector_rn', 'camera.dark_current')

        snr2 = -(_snr**2)
        fstar = self._fsource
        fsky = self.camera._fsky(verbose=self.verbose)
        Npix = self.camera._sn_box(self.verbose)
        thermal = self.camera.c_thermal(verbose=self.verbose)

        dark_rate = _dark_current[0] * u.Unit(_dark_current[1]) #<<-'electron / (pix s)'
        rn = _detector_rn[0] * u.Unit(_detector_rn[1])

        QE = _total_qe[0] * u.Unit(_total_qe[1])
        a = (QE * fstar)**2
        b = snr2 * (QE * (fstar + fsky) + thermal + dark_rate * Npix)
        c = snr2 * rn**2 * Npix * _nexp
        texp = ((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).to(u.s)

        self._exptime = texp

        return True

    def _update_magnitude(self):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """
        self.camera._print_initcon(self.verbose)

        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')
        (f0, c_ap, D, dlam) = self.recover('camera.ab_zeropoint',
                                           'camera.ap_corr',
                                           'telescope.effective_aperture',
                                           'camera.derived_bandpass')
        (_total_qe, _detector_rn, _dark_current) = self.recover('camera.total_qe',
                                    'camera.detector_rn',
                                    'camera.dark_current')

        exptime = (_exptime[0] * u.Unit(_exptime[1])).to(u.s)

        D = D.to(u.cm)
        fsky = self.camera._fsky(verbose=self.verbose)

        Npix = self.camera._sn_box(self.verbose)
        c_t = self.camera.c_thermal(verbose=self.verbose)

        QE = _total_qe[0] * u.Unit(_total_qe[1])
        rn = _detector_rn[0] * u.Unit(_detector_rn[1])
        dark_rate = _dark_current[0] * u.Unit(_dark_current[1]) #<<-'electron / (pix s)'

        snr2 = -(_snr ** 2)

        a0 = (QE * exptime)**2
        b0 = snr2 * QE * exptime
        c0 = snr2 * ((QE * fsky + c_t + Npix * dark_rate) * exptime + (rn**2 * Npix * _nexp))
        k = (-b0 + np.sqrt(b0**2 - 4. * a0 * c0)) / (2. * a0)
        flux = (4. * k) / (f0 * c_ap[0] * np.pi * D**2 * (dlam*u.nm))

        self._magnitude = -2.5 * np.log10(np.array(flux)) * u.mag('AB')

        return True

    def _update_snr(self):
        """
        Calculate the SNR for the given exposure time and source SED.
        """

        self.camera._print_initcon(self.verbose)

        (_exptime, _nexp, n_bands) = self.recover('_exptime', 'n_exp',
                                                  'camera.n_bands')

        (_total_qe, _detector_rn, _dark_current) = self.recover('camera.total_qe',
                             'camera.detector_rn', 'camera.dark_current')


        number_of_exposures = np.full(n_bands, _nexp)
        desired_exp_time = (np.full(n_bands, _exptime[0]) * u.Unit(_exptime[1])).to(u.second)
        time_per_exposure = desired_exp_time / number_of_exposures

        QE = _total_qe[0] * u.Unit(_total_qe[1])

        signal_counts = QE * self._fsource * desired_exp_time
        shot_noise_in_signal = np.sqrt(signal_counts)

        sky_counts = QE * self.camera._fsky(verbose=self.verbose) * desired_exp_time
        shot_noise_in_sky = np.sqrt(sky_counts)

        sn_box = self.camera._sn_box(self.verbose) #<-- units should be "pix"

        rn = _detector_rn[0] * u.Unit(_detector_rn[1])
        read_counts = rn**2 * sn_box * number_of_exposures

        dark_rate = _dark_current[0] * u.Unit(_dark_current[1])
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
            self._update_snr()
        if self.unknown == "exptime":
            self._update_exptime()

    def _update_snr(self):
        """
        Calculate the SNR based on the current SED and spectrograph parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('aperture'))
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
        if self.source.sed.fluxunits.name == "abmag":
            funit = u.ABmag
        elif self.source.sed.fluxunits.name == "photlam":
            funit = u.ph / u.s / u.cm**2 / u.AA
        elif self.source.sed.fluxunits.name == "flam":
            funit = u.erg / u.s / u.cm**2 / u.AA
        else:
            funit = u.Unit(self.source.sed.fluxunits.name)
        wave = _wave.to(u.AA)

        swave = (self.source.sed.wave * u.Unit(self.source.sed.waveunits.name)).to(u.AA)

        sflux = (self.source.sed.flux * funit).to(u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(swave))
        wave = wave.to(swave.unit)

        delta_lambda = self.recover('spectrograph.delta_lambda').to(u.AA / u.pix)

        iflux = np.interp(wave, swave, sflux, left=0., right=0.)
        self._interp_flux = iflux
        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.ct
        #print('_update_snr phot_energy: ', phot_energy)

        scaled_aeff = aeff * (aper / (15 * u.m))**2
        source_counts = iflux / phot_energy * scaled_aeff * exptime * delta_lambda
        #print('_update_snr source_counts: ', source_counts)

        bg_counts = bef / phot_energy * scaled_aeff * exptime

        snr = source_counts / np.sqrt(source_counts + bg_counts)

        if self.verbose:
            print("SNR: {}".format(snr))

        self._snr = snr

    def _update_exptime(self):
        """
        Calculate the exptime based on the current SED and spectrograph parameters.
        """

        if self.verbose:
            msg1 = "Creating exposure for {} ({})".format(self.telescope.name,
                                                           self.telescope.recover('aperture'))
            msg2 = " with {} in mode {}".format(self.spectrograph.name, self.spectrograph.mode)
            print(msg1 + msg2)

        _snr_goal, _exptime = self.recover('_snr_goal', '_exptime')
        _wave, aeff, bef, aper, R, wrange = self.recover('spectrograph.wave',
                                                         'spectrograph.aeff',
                                                         'spectrograph.bef',
                                                         'telescope.effective_aperture',
                                                         'spectrograph.R',
                                                         'spectrograph.wrange')

        if self.verbose:
            print("The requested SNR is {}\n".format(_snr_goal))

        if self.source.sed.fluxunits.name == 'abmag'  == "abmag":
            funit = u.ABmag
        elif self.source.sed.fluxunits.name  == "photlam":
            funit = u.ph / u.s / u.cm**2 / u.AA
        else:
            funit = u.Unit(self.source.sed.fluxunits.name)

        wave = _wave.to(u.AA)

        swave = (self.source.sed.wave * u.Unit(self.source.sed.waveunits.name)).to(u.AA)

        sflux = (self.source.sed.flux * funit).to(u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(swave))

        wave = wave.to(swave.unit)

        delta_lambda = self.recover('spectrograph.delta_lambda').to(u.AA / u.pix)

        iflux = np.interp(wave, swave, sflux, left=0., right=0.)

        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / wave.to(u.cm) / u.ct

        scaled_aeff = aeff * (aper / (15 * u.m))**2

        if (self.verbose):
            print('sflux  = ', sflux, '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
            print('wave = ', wave, '\n') #<---- 20,600 element array of wavelengths tied to Spectrograph object (not Exposure)
            print('delta_lambda = ', delta_lambda, '\n') #<--- this has the correct units, "Angstrom/pix"
            print('iflux = ', iflux, '\n') #<--- this has the correct units, "erg / (Angstrom s cm2)"
                                 #<--- becuase the units are carried through the interpolation
            print('bef = ', bef)  #<--- this has the correct units, "erg / (pix s cm2)"
            print('photE = ', phot_energy, '\n') #<--- this has the correct units, "erg / ct"
            print('aeff = ', aeff) #<--- this has the correct units, "cm2"
            print('aper = ', aper)#<--- this has the correct units, "m"
            print('scaled_aeff = ', scaled_aeff, '\n') #<--- this has the correct units, "cm2"
            print('SNR^2 :', (_snr_goal)**2)

        t_exp = (_snr_goal)**2 * (iflux / phot_energy * scaled_aeff * delta_lambda + bef / phot_energy * scaled_aeff) / ((iflux/phot_energy)**2 * scaled_aeff**2 * delta_lambda**2)

        if self.verbose:
            print("Exptime: {}".format(t_exp))

        self._exptime = t_exp

        return True #completed successfully

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
        print('doesnt exist yet pull it from camera class')

        return False #completed successfully

    def _update_magnitude(self):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """

        print('doesnt exist yet pull it from camera class')

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
