
def uvspec_snr(telescope, band, template, fuvmag, exptime, silent=False):
    ''' Run a basic SNR calculation that takes in a telescope,
        spectral template, normalization magnitude, and exposure
        time to compute SNR. For converting magnitude, template,
	      and SNR to a desired exposure time, use uvspec_exptime.py

        usage:
	      wave, snr, uvi = uvspec_snr(telescope, band, template, uvmag, exptime)

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string.
             EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis
             EAC2 = 6 m diameter off-axis
             EAC3 = 8 m diameter on-axis

           2-band = your choice of UVI grating, a string:
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
    import numpy as np
    import astropy.units as u

    # create the basic objects 
    tel = Telescope()
    tel.set_from_hwome(telescope)
    suitable_instruments, suitable_bands = tel.find_instrument_with("disperser")
    instrument = None
    for test_band in suitable_bands:
        if band in test_band:
            instrument = suitable_bands[test_band]
            break
    if instrument is None:
        raise ValueError(f"Could not find an instrument with {band}")
    print(f"Using Instrument {instrument} with band {test_band}")

    source = Source()
    redshift = 0.0
    extinction = 0.0
    source.set_sed(template, fuvmag, redshift, extinction, bandpass="galex,fuv")

    tel.verbose = True
    inst = tel.instruments[instrument]

    uvi_exp = SourceSpectrographicExposure()
    uvi_exp.source = source
    uvi_exp.verbose = not silent

    inst.add_exposure(uvi_exp)
    inst.band = test_band # doing it this way is a little more forgiving as an API

    if (silent):
        uvi_exp.verbose = False
        tel.verbose = False
        inst.verbose = False
        print("We have set verbose = False")

    if not silent:
        print("Current SED template: {}".format(uvi_exp.source.name))
        print("Current grating mode: {}".format(inst.band))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp.exptime = exptime * u.hr

    uvi_exp.unknown = "snr"

    uvi_snr = uvi_exp.recover('snr')

    wave, snr =  uvi_exp.wave, uvi_exp.snr


    return wave, snr, inst



def uvspec_exptime(telescope, band, template, fuvmag, snr, silent=False):

    ''' 
    Run a basic SNR calculation that takes in a telescope, spectral template,
    normalization magnitude, and SNR goal to compute exposure time. For converting
    magnitude, template, and exptime to SNR, use uvspec_snr.py

      usage:
	      wave, exptime, uvi = uvspec_exptime(telescope, band, template, uvmag, snr)

        positional arguments:

          1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string.
            EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
            EAC2 = 6 m diameter off-axis 
            EAC3 = 8 m diameter on-axis

          2-band = your choice of UVI grating, a string:
          ['G120M', 'G150M', 'G180M', 'G155L', 'G145LL', 'G300M']

          3-template = your choice of spectral template:
            ['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts',
                      'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']

          4-fuvmag = FUV magnitude to normalize the template spectrum, a float.

          5-snr = desired SNR, per pixel

      outputs are two arrays of floats for wavelength and exptime and the Spectrograph
      object in case it is needed by other code.
    '''

    from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure
    import astropy.units as u

    # create the basic objects
    tel = Telescope()
    tel.set_from_hwome(telescope)
    suitable_instruments, suitable_bands = tel.find_instrument_with("disperser")
    instrument = None
    for test_band in suitable_bands:
        if band in test_band:
            instrument = suitable_bands[test_band]
            break
    if instrument is None:
        raise ValueError(f"Could not find an instrument with {band}")
    print(f"Using Instrument {instrument} with band {test_band}")


    source = Source()
    redshift = 0.0
    extinction = 0.0
    source.set_sed(template, fuvmag, redshift, extinction, bandpass="galex,fuv")

    uvi_exp = SourceSpectrographicExposure() 
    uvi_exp.source = source
    uvi_exp.verbose = not silent

    inst = tel.instruments[instrument]

    uvi_exp = SourceSpectrographicExposure() 
    uvi_exp.source = source
    uvi_exp.verbose = not silent

    inst.add_exposure(uvi_exp)
    inst.band = test_band  # doing it this way is a little more forgiving as an API

    if not silent:
        print("Current SED template: {}".format(template))
        print("Current grating mode: {}".format(inst.band))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp.snr = snr# * u.ct**0.5 / u.pix**0.5
    uvi_exp.unknown = 'exptime' #< --- this triggers the _update_exptime function in the SpectrographicExposure exposure object

    uvi_exptime = uvi_exp.recover('exptime')

    wave, exptime =  uvi_exp.wave, uvi_exp.exptime


    return wave, exptime, inst
