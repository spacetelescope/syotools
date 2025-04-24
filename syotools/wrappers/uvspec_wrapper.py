
def uvspec_snr(telescope, mode, template, fuvmag, exptime, silent=False):
    ''' Run a basic SNR calculation that takes in a telescope,
        spectral template, normalization magnitude, and exposure
        time to compute SNR. For converting magnitude, template,
	      and SNR to a desired exposure time, use uvspec_exptime.py

        usage:
	      wave, snr, uvi = uvspec_snr(telescope, mode, template, uvmag, exptime)

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string.
             EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis
             EAC2 = 6 m diameter off-axis
             EAC3 = 8 m diameter on-axis

           2-mode = your choice of UVI grating, a string:
		        ['G120M', 'G150M', 'G180M', 'G155L', 'G145LL', 'G300M']

           3-template = your choice of spectral template:
		          ['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts',
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']

           4-fuvmag = FUV magnitude to normalize the template spectrum, a float.

	   5-exptime = desired exposure time in hours, a float

        outputs are two arrays of floats for wavelength and snr and the Spectrograph
		    object in case it is needed by other code.
    '''

    from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure
    import numpy as np, astropy.units as u

    # create the basic objects 
    uvi, tel = Spectrograph(), Telescope() 
    tel.set_from_sei(telescope)
    tel.add_spectrograph(uvi)
    uvi.mode = mode

    source = Source()
    redshift = 0.0
    extinction = 0.0
    source.set_sed(template, fuvmag, redshift, extinction, bandpass="galex,fuv")

    uvi_exp = SourceSpectrographicExposure()
    uvi_exp.source = source
    uvi_exp.verbose = not silent
    uvi.add_exposure(uvi_exp)

    uvi_exp.exptime = [[exptime, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], 'hr']

    tel.verbose = True
    if (silent):
       uvi_exp.verbose = False
       tel.verbose = False
       uvi.verbose = False
       print("We have set verbose = False")

    if not silent:
        print("Current SED template: {}".format(uvi_exp.source.name))
        print("Current grating mode: {}".format(uvi.descriptions[uvi.mode]))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp.enable()
    uvi_snr = uvi_exp.recover('snr')

    wave, snr =  uvi.wave, uvi_exp.snr

    return wave, snr, uvi



def uvspec_exptime(telescope, mode, template, fuvmag, snr_goal, silent=False):

    ''' Run a basic SNR calculation that takes in a telescope,
      spectral template, normalization magnitude, and SNR goal
      to compute exposure time. For converting magnitude, template,
      and exptime to SNR, use uvspec_snr.py

        usage:
	       wave, exptime, uvi = uvspec_exptime(telescope, mode, template, uvmag, snr_goal)

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string.
             EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis
             EAC2 = 6 m diameter off-axis
             EAC3 = 8 m diameter on-axis

           2-mode = your choice of UVI grating, a string:
		        ['G120M', 'G150M', 'G180M', 'G155L', 'G145LL', 'G300M']

           3-template = your choice of spectral template:
	          	['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts',
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']

           4-fuvmag = FUV magnitude to normalize the template spectrum, a float.

	         5-snr_goal = desired SNR, per pixel

         outputs are two arrays of floats for wavelength and exptime and the Spectrograph
		     object in case it is needed by other code.
       '''

    from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure
    import astropy.units as u

    # create the basic objects
    uvi, tel = Spectrograph(), Telescope()
    tel.set_from_json(telescope)
    tel.add_spectrograph(uvi)
    uvi.mode = mode

    source = Source()
    redshift = 0.0
    extinction = 0.0
    source.set_sed(template, fuvmag, redshift, extinction, bandpass="galex,fuv")

    uvi_exp = SourceSpectrographicExposure()
    uvi_exp.source = source
    uvi_exp.verbose = not silent
    uvi.add_exposure(uvi_exp)

    if not silent:
        print("Current SED template: {}".format(source.sed.name))
        print("Current grating mode: {}".format(uvi.descriptions[uvi.mode]))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp._snr_goal= snr_goal * (u.ct)**0.5 / (u.pix)**0.5

    snr = uvi_exp.recover('exptime')
    uvi_exp.unknown = 'exptime' #< --- this triggers the _update_exptime function in the SpectrographicExposure exposure object

    uvi_exptime = uvi_exp.recover('exptime')

    wave, exptime =  uvi.wave, uvi_exp.exptime

    return wave, exptime, uvi
