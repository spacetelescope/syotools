def camera_snr(telescope, template, mag, exptime, silent=False): 

    ''' Run a basic SNR calculation that takes in a telescope, 
      spectral template, normalization magnitude, and exptime   
      to compute SNR. For converting magnitude, template, 
      and SNR to exptime, use camera_exptime.py 
      
        usage: 
           snr, hri = camera_snr(telescope, template, mag, snr_goal) 

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
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

           3-mag = V magnitude to normalize the template spectrum, a float. 

           4-exptime =  desired exptime in hours 

           outputs are arrays with the SNR in each band for FUV, NUV, U, B, V, R, I, J, H, K 
           and the camera object "hri" 
       '''

    
    from syotools.models import Telescope, Camera
    from syotools.utils.jsonunit import str_jsunit
    from syotools.spectra import SpectralLibrary
    import astropy.units as u
    import numpy as np 

    # create the basic objects 
    hri, tel = Camera(), Telescope() 
    tel.set_from_sei(telescope)
    tel.add_camera(hri)
    hri_exp = hri.create_exposure()

    hri_exp.sed_id = template
    hri_exp.renorm_sed(mag * u.ABmag, bandpass='v')
    exptime_list = []
    for i in hri.bandnames: exptime_list.append(exptime) 
    hri_exp.exptime[1]['value'] = exptime_list 

    #Print the current template & mode
    if not silent: 
        print("Current SED template: {}".format(hri_exp.sed_id)) 
        print("Current exposure time: {} hours\n".format(hri_exp.exptime[1]['value'][0])) 
    
    hri_exp.enable()
    snr = hri_exp.recover('snr')
    hri_sed, hri_snr = hri_exp.recover('sed', 'snr')  
    hri_snr = hri_exp.snr[1]['value']

    if not silent: 
        for bb, ss in zip(hri.bandnames, hri_snr): print("{}, SNR = {}".format(bb, ss)) 

    return snr, hri_snr



def camera_exptime(telescope, template, mag, snr_goal, silent=False): 

    ''' Run a basic SNR calculation that takes in a telescope, 
      spectral template, normalization magnitude, and SNR goal  
      to compute exptime. For converting magnitude, template, 
      and exptime to SNR, use camera_snr.py 
      
        usage: 
          exptime, hri = camera_exptime(telescope, template, mag, snr_goal) 

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
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

           3-mag = V magnitude to normalize the template spectrum, a float.

           4-snr_goal = desired SNR, per pixel, for each band 

           outputs are arrays with the SNR in each band for FUV, NUV, U, B, V, R, I, J, H, K 
       '''

    from syotools.models import Telescope, Camera
    from syotools.utils.jsonunit import str_jsunit
    from syotools.spectra import SpectralLibrary
    import astropy.units as u
    import numpy as np

    # create the basic objects 
    hri, tel = Camera(), Telescope()
    tel.set_from_sei(telescope)
    tel.add_camera(hri)
    hri_exp = hri.create_exposure()

    hri_exp.sed_id = template
    hri_exp.renorm_sed(mag * u.ABmag, bandpass='v')
    snr_list = []
    for i in hri.bandnames: snr_list.append(snr_goal) 
    hri_exp._snr[1]['value'] = snr_list 

    #Print the current template & SNR goal
    if not silent: 
        print("Current SED template: {}".format(hri_exp.sed_id))
        print("Current SNR goals: {}".format(hri_exp._snr[1]["value"])) 
        
    hri_exp.unknown = 'exptime'
    hri_exptime = hri_exp.recover('exptime')

    if not silent: 
        for bb, ee in zip(hri.bandnames, hri_exptime): print("{}, exptime = {}".format(bb, ee))
    
    return hri_exptime, hri
