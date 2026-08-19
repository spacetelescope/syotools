def camera_snr(telescope, template, magnitude, exptime, silent=False): 
	''' 
	Run a basic SNR calculation that takes in a telescope, spectral template,
	normalization magnitude, and exptime to compute SNR.
	
	usage: snr, hri = camera_snr(telescope, template, mag, exptime) 
	
	positional arguments:
	
	1 - telescope = 'EAC5'. This argument is a string. 
	    EAC5 = 10 m outer diameter mirror
	
	2 - spectral template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 'Seyfert
	    1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 'G191B2B (WD)',
	    'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 'Starburst, E(B-V) = 0.6', 'B5V
	    Star', 'M2V Star', 'Elliptical Galaxy', 'Sbc Galaxy', 'Starburst Galaxy', 'NGC
	    1068', 'Galaxy with f_esc, HI=1, HeI=1', 'Galaxy with f_esc, HI=0.001, HeI=1',
	    'Blackbody5000', 'Blackbody100000' 

	3 - mag = V magnitude to normalize the template spectrum

	4 - exptime = desired exptime in hours 
	
	outputs are dicts of the snrs, and of the instrument objects, keyed by the 
    filter name
    '''

	from syotools.models import Telescope, Source, SourcePhotometricExposure
	import numpy as np
	import astropy.units as u 
      
	tel = Telescope()   # create a Telescope, Camera, and Exposure 
	tel.set_from_hwome(telescope)
	suitable_instruments, suitable_filters = tel.find_instrument_with("filter")
	
	source = Source() 
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, magnitude, redshift, extinction, bandpass="johnson,v")

	results = {}
	out_instrument = {}

	for instrument in suitable_instruments:
		inst = tel.instruments[instrument]

		exp = SourcePhotometricExposure() 
		exp.source = source

		inst.add_exposure(exp)
		exp.exptime = [exptime] * u.hr
		exp.unknown = 'snr'

		if not silent: 
			print('------ Computing SNR as the Unknown -------') 
			for bb, ss in zip(inst.bands, exp.snr): print("{}, SNR = {}".format(bb, ss)) 
	            
		for bb,ee in zip(inst.bands, exp.snr): 
			results[bb] = ee
			out_instrument[bb] = inst

	return results, out_instrument


def camera_exptime(telescope, template, magnitude, snr, silent=False): 
	''' 
	Run a basic SNR calculation that takes in a telescope, spectral template,
	normalization magnitude, and SNR goal to compute exposure time. 
	
	usage: exptime, inst = camera_exptime(telescope, template, mag, snr) 
	
	positional arguments:
	
	1 - telescope = 'EAC5'. This argument is a string. 
	    EAC5 = 10 m outer diameter mirror
	
	2 - template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 'Seyfert
	    1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 'G191B2B (WD)',
	    'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 'Starburst, E(B-V) = 0.6', 'B5V
	    Star', 'M2V Star', 'Elliptical Galaxy', 'Sbc Galaxy', 'Starburst Galaxy', 'NGC
	    1068', 'Galaxy with f_esc, HI=1, HeI=1', 'Galaxy with f_esc, HI=0.001, HeI=1',
	    'Blackbody5000', 'Blackbody100000' 
	
	3 - magnitude = V magnitude to normalize the template spectrum, a float.
	
	4 - snr = desired SNR, per pixel, for each band
	
	outputs are dicts of the snrs, and of the instrument objects, keyed by the 
    filter name
	'''

	from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
	import numpy as np
	import astropy.units as u 
    
	# create a Telescope, Camera, and Exposure 
	tel = Telescope()
	tel.set_from_hwome(telescope)
	suitable_instruments, suitable_filters = tel.find_instrument_with("filter")
	
	source = Source()
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, magnitude, redshift, extinction, bandpass="johnson,v")

	results = {}
	out_instrument = {}

	for instrument in suitable_instruments:
		inst = tel.instruments[instrument]
		exp = SourcePhotometricExposure()
		exp.source = source

		inst.add_exposure(exp)

		exp.snr = snr #* u.Unit('electron(1/2)')

		exp.unknown = 'exptime'
		if not silent: 
			print('-- Computing Exptime as the Unknown --') 
			for bb, ee in zip(inst.bands, exp.exptime): print("{}, exptime = {}".format(bb, ee))

		for bb,ee in zip(inst.bands, exp.exptime): 
			results[bb] = ee
			out_instrument[bb] = inst

	return results, out_instrument

def camera_magnitude(telescope, template, snr, exptime, silent=False): 
	''' 
	Run a basic magnitude calculation that takes in a telescope, spectral template,
	snr, and exptime to compute magnitude.
	
	usage: mag, instrument = camera_magnitude(telescope, template, snr, exptime) 
	
	positional arguments:
	
	1 - telescope = 'EAC5'. This argument is a string. 
	    EAC5 = 10 m outer diameter mirror
	
	2 - spectral template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 'Seyfert
	    1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 'G191B2B (WD)',
	    'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 'Starburst, E(B-V) = 0.6', 'B5V
	    Star', 'M2V Star', 'Elliptical Galaxy', 'Sbc Galaxy', 'Starburst Galaxy', 'NGC
	    1068', 'Galaxy with f_esc, HI=1, HeI=1', 'Galaxy with f_esc, HI=0.001, HeI=1',
	    'Blackbody5000', 'Blackbody100000' 

	3 - snr = desired SNR, per pixel, for each band

	4 - exptime = desired exptime in hours 
	
	outputs are dicts of the magnitudes, and of the instrument objects, keyed by the 
    filter name
    '''

	from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
	import numpy as np, astropy.units as u

	tel = Telescope()
	# create a Telescope, Camera, and Exposure 
	tel.set_from_hwome(telescope)
	suitable_instruments, suitable_filters = tel.find_instrument_with("filter")
	
	source = Source() 
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, 30., redshift, extinction)   

	results = {}
	out_instrument = {}

	for instrument in suitable_instruments:
		inst = tel.instruments[instrument]
		exp = SourcePhotometricExposure()
		exp.source = source

		inst.add_exposure(exp)
		exp.exptime = [exptime] * u.hr
		exp.snr = [snr] * u.dimensionless_unscaled
		exp.unknown = 'magnitude' 

		if not silent: 
			print('--- Computing Magnitude as the Unknown ---') 
			for bb, mm in zip(inst.bands, exp.magnitude): print("{}, Mag = {}".format(bb, mm)) 
	
		for bb,ee in zip(inst.bands, exp.magnitude): 
			results[bb] = ee
			out_instrument[bb] = inst

	return results, out_instrument
