
def uvspec_snr(telescope, band, template, fuvmag, exptime, silent=False):
    ''' 
    Run a basic SNR calculation that takes in a telescope,
    spectral template, normalization magnitude, and exposure
    time to compute SNR. For converting magnitude, template,
        and SNR to a desired exposure time, use uvspec_exptime.py

    usage:
        wave, snr, inst = uvspec_snr(telescope, band, template, uvmag, exptime)

        positional arguments:

	1 - telescope = 'EAC5'. This argument is a string. 
	    EAC5 = 10 m outer diameter mirror

    2 - band = your choice of grating, a string:
            ['HRI_S_NIR.HRI_Grism1a', 'HRI_S_NIR.HRI_Grism1b', 
             'HRI_S_UVIS.HRI_Grism_2a', 'HRI_S_UVIS.HRI_Grism_2b', 
             'FUV_MOS.G120M', 'FUV_MOS.G150M', 'FUV_MOS.G180M', 
             'FUV_MOS_L.G145LL', 'FUV_MOS_L.G155L', 'FUV_MOS_L.G165LL', 
             'NUV_MOS.G300M', 'NUV_MOS.G700L']

	3 - template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 'Seyfert
	    1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 'G191B2B (WD)',
	    'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 'Starburst, E(B-V) = 0.6', 'B5V
	    Star', 'M2V Star', 'Elliptical Galaxy', 'Sbc Galaxy', 'Starburst Galaxy', 'NGC
	    1068', 'Galaxy with f_esc, HI=1, HeI=1', 'Galaxy with f_esc, HI=0.001, HeI=1',
	    'Blackbody5000', 'Blackbody100000' 

    4 - fuvmag = FUV magnitude to normalize the template spectrum, a float.

    5 - exptime = desired exposure time in hours, a float

    outputs are two arrays of floats for wavelength and snr, and the Spectrograph
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
    # this code demonstrates how to find a band with a partial name
    for test_band in suitable_bands:
        if band in test_band:
            instrument = suitable_bands[test_band]
            break
    if instrument is None:
        raise ValueError(f"Could not find an instrument with {band}")

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
        print(f"Using Instrument {instrument} with band {test_band}")
        print("Current SED template: {}".format(uvi_exp.source.name))
        print("Current grating mode: {}".format(inst.band))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp.exptime = exptime * u.hr

    uvi_exp.unknown = "snr"

    uvi_snr = uvi_exp.recover('snr')

    wave, snr =  uvi_exp.wave, uvi_exp.snr[0]


    return wave, snr, inst



def uvspec_exptime(telescope, band, template, fuvmag, snr, silent=False):

    ''' 
    Run a basic SNR calculation that takes in a telescope, spectral template,
    normalization magnitude, and SNR goal to compute exposure time. For converting
    magnitude, template, and exptime to SNR, use uvspec_snr.py

      usage:
	      wave, exptime, inst = uvspec_exptime(telescope, band, template, uvmag, snr)

        positional arguments:

	1 - telescope = 'EAC5'. This argument is a string. 
	    EAC5 = 10 m outer diameter mirror

    2 - band = your choice of grating, a string:
            ['HRI_S_NIR.HRI_Grism1a', 'HRI_S_NIR.HRI_Grism1b', 
             'HRI_S_UVIS.HRI_Grism_2a', 'HRI_S_UVIS.HRI_Grism_2b', 
             'FUV_MOS.G120M', 'FUV_MOS.G150M', 'FUV_MOS.G180M', 
             'FUV_MOS_L.G145LL', 'FUV_MOS_L.G155L', 'FUV_MOS_L.G165LL', 
             'NUV_MOS.G300M', 'NUV_MOS.G700L']

	3 - template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 'Seyfert
	    1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 'G191B2B (WD)',
	    'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 'Starburst, E(B-V) = 0.6', 'B5V
	    Star', 'M2V Star', 'Elliptical Galaxy', 'Sbc Galaxy', 'Starburst Galaxy', 'NGC
	    1068', 'Galaxy with f_esc, HI=1, HeI=1', 'Galaxy with f_esc, HI=0.001, HeI=1',
	    'Blackbody5000', 'Blackbody100000' 

    4 - fuvmag = FUV magnitude to normalize the template spectrum, a float.

    5 - snr = desired SNR, per pixel

    outputs are two arrays of floats for wavelength and snr, and the Spectrograph
        object in case it is needed by other code.
    '''
    from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure
    import astropy.units as u

    # create the basic objects
    tel = Telescope()
    tel.set_from_hwome(telescope)
    suitable_instruments, suitable_bands = tel.find_instrument_with("disperser")
    instrument = None
    # this code demonstrates how to find a band with a partial name
    for test_band in suitable_bands:
        if band in test_band:
            instrument = suitable_bands[test_band]
            break
    if instrument is None:
        raise ValueError(f"Could not find an instrument with {band}")


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
        print(f"Using Instrument {instrument} with band {test_band}")
        print("Current SED template: {}".format(template))
        print("Current grating mode: {}".format(inst.band))
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime))

    uvi_exp.snr = snr# * u.ct**0.5 / u.pix**0.5
    uvi_exp.unknown = 'exptime' #< --- this triggers the _update_exptime function in the SpectrographicExposure exposure object

    uvi_exptime = uvi_exp.recover('exptime')

    wave, exptime =  uvi_exp.wave, uvi_exp.exptime[0]


    return wave, exptime, inst
