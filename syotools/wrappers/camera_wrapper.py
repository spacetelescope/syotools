def camera_snr(telescope, template, magnitude, exptime, silent=False): 
	''' Run a basic SNR calculation that takes in a telescope, 
	spectral template, normalization magnitude, and exptime   
	to compute SNR.
	
	usage: 
	snr, hri = camera_snr(telescope, template, mag, snr_goal) 
	
	positional arguments:
	
	telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 

	    EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
	    EAC2 = 6 m diameter off-axis 
	    EAC3 = 8 m diameter on-axis 
	
	spectral template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 
	    'Seyfert 1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 
	    'G191B2B (WD)', 'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 
	    'Starburst, E(B-V) = 0.6', 'B5V Star', 'M2V Sta', 'Elliptical Galaxy', 
	    'Sbc Galaxy', 'Starburst Galaxy', 'NGC 1068', 'Galaxy with f_esc, HI=1, HeI=1', 
	    'Galaxy with f_esc, HI=0.001, HeI=1', 'Blackbody5000', 'Blackbody100000' 

	mag = V magnitude to normalize the template spectrum

	exptime = desired exptime in hours 
	
	outputs are arrays with the SNR in each band for FUV, NUV, U, B, V, R, I, J, H, K 
	and the camera object "hri" 
       '''

	from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
	import numpy as np, astropy.units as u 
      
	tel, hri = Telescope(), Camera()   # create a Telescope, Camera, and Exposure 
	tel.set_from_sei(telescope)
	hri.set_from_sei('HRI')
	
	source = Source() 
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, magnitude, redshift, extinction, bandpass="johnson,v")   
	
	exp = SourcePhotometricExposure() 
	exp.source = source
	    
	exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']
	exp.unknown = 'snr'
	tel.add_camera(hri)
	hri.add_exposure(exp)

	if not silent: 
		print('------ Computing SNR as the Unknown -------') 
		for bb, ss in zip(hri.bandnames, exp.snr): print("{}, SNR = {}".format(bb, ss)) 
	            
	return exp.snr, hri 


def camera_exptime(telescope, template, magnitude, snr_goal, silent=False): 
	''' Run a basic SNR calculation that takes in a telescope, 
	spectral template, normalization magnitude, and SNR goal  
	to compute exposure time. 
	
	usage: 
	exptime, hri = camera_exptime(telescope, template, mag, snr_goal) 
	
	positional arguments:
	
	telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
	    EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
	    EAC2 = 6 m diameter off-axis 
	    EAC3 = 8 m diameter on-axis 
	
	2-template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 
	    'Seyfert 1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 
	    'G191B2B (WD)', 'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 
	    'Starburst, E(B-V) = 0.6', 'B5V Star', 'M2V Sta', 'Elliptical Galaxy', 
	    'Sbc Galaxy', 'Starburst Galaxy', 'NGC 1068', 'Galaxy with f_esc, HI=1, HeI=1', 
	    'Galaxy with f_esc, HI=0.001, HeI=1', 'Blackbody5000', 'Blackbody100000' 
	
	magnitude = V magnitude to normalize the template spectrum, a float.
	
	snr_goal = desired SNR, per pixel, for each band 
	
	outputs	are arrays with the exptime in each band for FUV, NUV, U, B, V, R, I, J, H, K 
	'''

	from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
	import numpy as np, astropy.units as u 
    
	# create a Telescope, Camera, and Exposure 
	tel, hri = Telescope(), Camera()
	tel.set_from_sei(telescope)
	hri.set_from_sei('HRI')
	
	source = Source() 
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, magnitude, redshift, extinction, bandpass="johnson,v")   

	exp = SourcePhotometricExposure() 
	exp.source = source
	
	exp._snr = [snr_goal] * u.Unit('electron(1/2)')  
	exp.unknown = 'exptime' 
	tel.add_camera(hri)
	hri.add_exposure(exp)
	
	if not silent: 
		print('-- Computing Exptime as the Unknown --') 
		for bb, ee in zip(hri.bandnames, exp.exptime): print("{}, SNR = {}".format(bb, ee)) 

	return exp.exptime, hri 

def camera_magnitude(telescope, template, snr, exptime, silent=False): 
	''' Run a basic SNR calculation that takes in a telescope, 
	spectral template, SNR goal, and exposure time and computes 
	the limiting magnitude. 
	
	usage: 
	exptime, hri = camera_magnitude(telescope, template, mag, snr_goal) 
	
	positional arguments:
	
	telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
	    EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
    	    EAC2 = 6 m diameter off-axis 
	    EAC3 = 8 m diameter on-axis 

	spectral template = your choice of spectral template: 
	    'Classical T Tauri', 'M1 Dwarf', 'G Dwarf', '10 Myr Starburst', 'QSO', 
	    'Seyfert 1', 'Seyfert 2', 'Liner', 'O5V Star', 'G2V Star', 'Orion Nebula', 
	    'G191B2B (WD)', 'GD71 (WD)', 'GD153 (WD)', 'Starburst, No Dust', 
	    'Starburst, E(B-V) = 0.6', 'B5V Star', 'M2V Sta', 'Elliptical Galaxy', 
	    'Sbc Galaxy', 'Starburst Galaxy', 'NGC 1068', 'Galaxy with f_esc, HI=1, HeI=1', 
	    'Galaxy with f_esc, HI=0.001, HeI=1', 'Blackbody5000', 'Blackbody100000' 

	snr_goal = desired SNR for each band 

	exptime = exposure time per band in hours

	outputs are arrays with the limiting in each band for FUV, NUV, U, B, V, R, I, J, H, K 
	'''

	from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
	import numpy as np, astropy.units as u
	
	# create a Telescope, Camera, and Exposure 
	tel, hri = Telescope(), Camera()
	tel.set_from_sei(telescope)
	hri.set_from_sei('HRI')
	
	source = Source() 
	redshift = 0. # changes to these are not implemented yet 
	extinction = 0. 
	
	source.set_sed(template, 30., redshift, extinction)   
	        
	exp = SourcePhotometricExposure() 
	exp.source = source
	
	exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']
	exp._snr = [snr] * u.Unit('electron(1/2)')  
	    
	exp.unknown = 'magnitude' 
	tel.add_camera(hri)
	hri.add_exposure(exp)
	
	if not silent: 
		print('--- Computing Magnitude as the Unknown ---') 
		for bb, mm in zip(hri.bandnames, exp.magnitude): print("{}, SNR = {}".format(bb, mm)) 
	
	return exp.magnitude, hri 
