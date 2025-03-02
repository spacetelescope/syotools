;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;              ATLAST noise model - Tyler D. Robinson
;  inputs:
;    Ahr  - hi-res planetary albedo spectrum
;  lamhr  - wavelength grid for Ahr (um)
;  alpha  - phase angle (deg)
;    Phi  - phase function evaluated at alpha
;     Rp  - planetary radius (Rearth)
;   Teff  - stellar effective temperature (K)
;     Rs  - stellar radius (Rsun)
;      r  - orbital separation (au)
;      d  - distance to system (pc)
;    Nez  - number of exozodis
;  lammin - minimum wavelength (um)
;  lammax - maximum wavelength (um)
;     Res - spectral resolution (lambda/Dlambda)
;    diam - telescope diameter (m)
;    Tput - system throughput
;      C  - raw contrast
;     IWA - inner working angle (lambda/D)
;     OWA - outer working angle (lambda/D; unless /FIXOWA)
;    Tsys - observing system temperature for thermal background (K)
;    emis - observing system emissivity
;
;  outputs:
;    lam  - low-res wavelength grid (um)
;   dlam  - spectral bin widths (um)
;      A  - planetary albedo spectrum at low-res
;      q  - quantum efficiency
; Cratio  - planet-star contrast (flux) ratio
;     cp  - planet count rate (s**-1)
;    csp  - speckle count rate (s**-1)
;     cz  - zodiacal light count rate (s**-1) 
;    cez  - exozodiacal light count rate (s**-1)
;     cD  - dark current count rate (s**-1)
;     cR  - read noise count rate (s**-1)
;    cth  - internal thermal count rate (s**-1)
;  DtSNR  - integration time for SNR (hr)
;  wantsnr - SNR you want the integration time for
;
;  options:
;        FIX_OWA - set to fix OWA at OWA*lammin/D, as would occur 
;                  if lenslet array is limiting the OWA
;    COMPUTE_LAM - set to compute lo-res wavelength grid, otherwise 
;                  the grid input as variable 'lam' is used
;         SILENT - set to silence all warnings
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO ATLAST_NOISE, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                  lammin, lammax, Res, diam, Tput, C, IWA, OWA, Tsys, emis, $
                  lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, $
                  cth, DtSNR, wantsnr, whichplanet, FIX_OWA = fix_owa, COMPUTE_LAM = compute_lam, $
                  SILENT = silent

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  set key system params   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  De     = 1.d-4  ; dark current (s**-1)
  DNhpix = 3      ; horizontal pixel spread of IFS spectrum
  Re     = 0.1    ; read noise per pixel
  Dtmax  = 1.     ; maximum exposure time (hr)
  X      = 0.7    ; size of photometric aperture (lambda/D)
  q      = 0.9    ; quantum efficiency
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set astrophys parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  MzV  = 23.0 ; zodiacal light surface brightness (mag/arcsec**2)
  MezV = 22.0 ; exozodiacal light surface brightness (mag/arcsec**2)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  angular size of lenslet ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  theta = lammin/1.d6/diam/2.*(180/!DPI*3600.) ;assumes sampled at ~lambda/2D (arcsec)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set wavelength grid    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam) THEN BEGIN
    lam  = lammin ;in [um]
    Nlam = 1
    WHILE lam LT lammax DO BEGIN
      lam  = lam + lam/res
      Nlam = Nlam +1
    ENDWHILE 
    lam    = FLTARR(Nlam)
    lam[0] = lammin
    FOR j=1,Nlam-1 DO BEGIN
      lam[j] = lam[j-1] + lam[j-1]/res
    ENDFOR
  ENDIF
  Nlam = N_ELEMENTS(lam)
  dlam = FLTARR(Nlam) ;grid widths (um)
  FOR j=1,Nlam-2 DO BEGIN
    dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
  ENDFOR
  ;widths at edges are same as neighbor
  dlam[0] = dlam[1]
  dlam[Nlam-1] = dlam[Nlam-2]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      set throughput      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  T    = DBLARR(Nlam)
  T[*] = Tput
  sep  = r/d*SIN(alpha*!DPI/180.)*!DPI/180./3600.d ; separation in radians
  iIWA = WHERE( sep LT IWA*lam/diam/1.d6 )
  IF (iIWA[0] NE -1) THEN BEGIN
    T[iIWA] = 0. ;zero transmission for points inside IWA have no throughput
    IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: portions of spectrum inside IWA'
  ENDIF
  IF KEYWORD_SET(fix_owa) THEN BEGIN
    IF ( sep GT OWA*lammin/diam/1.d6 ) THEN BEGIN
      T[*] = 0. ;planet outside OWA, where there is no throughput
      IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: planet outside fixed OWA'
    ENDIF
  ENDIF ELSE BEGIN
    iOWA = WHERE( sep GT OWA*lam/diam/1.d6 )
    IF (iOWA[0] NE -1) THEN BEGIN
      T[iOWA] = 0. ;points outside OWA have no throughput
      IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: portions of spectrum outside OWA'
    ENDIF
  ENDELSE
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; degrade albedo spectrum  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam)  THEN A = DEGRADE_SPEC(Ahr,lamhr,lam,DLAM=dlam)
  IF ~KEYWORD_SET(compute_lam) THEN A = Ahr
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      compute fluxes      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Fs = Fstar(lam, Teff, Rs, r, /AU) ;stellar flux on planet
  Fp = Fplan(A, Phi, Fs, Rp, d)     ;planet flux at telescope
  Cratio = FpFs(A, Phi, Rp, r)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    compute count rates   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  cp     =  cplan(q, X, T, lam, dlam, Fp, diam)  ;planet count rate
  cz     =  czodi(q, X, T, lam, dlam, diam, MzV)  ;solar system zodi count rate
  cez    =  cezodi(q, X, T, lam, dlam, diam, r, Fstar(lam,Teff,Rs,1.,/AU), Nez, MezV)  ;exo-zodi count rate
  csp    =  cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam) ;speckle count rate
  cD     =  cdark(De, X, lam, diam, theta, DNhpix)  ;dark current count rate
  cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)  ;readnoise count rate
  cth    =  ctherm(q, X, lam, dlam, diam, Tsys, emis) ;internal thermal count rate
  cnoise =  cp + 2*(cz + cez + csp + cD + cR + cth) ; assumes background subtraction
  cb = (cz + cez + csp + cD + cR + cth)
  ctot = cp + cz + cez + csp + cD + cR + cth
  stop
  ;giada: where does the factor of 2 come from?
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  exposure time to SNR    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DtSNR = DBLARR(Nlam)
  DtSNR[*] = 0.d
  i     = WHERE(cp GT 0)
  IF (i[0] NE -1) THEN DtSNR[i] = (wantsnr^2.*cnoise[i])/cp[i]^2./3600. ; (hr)

  ;added by Giada:
  if whichplanet eq 'earth' then begin
  pt5 = closest(lam, 0.55) ;same wavelength chris stark used
  time = dtsnr(pt5)*3600.*1
  save, time, filename='earthtime.sav'
endif

  if whichplanet ne 'earth' then restore, 'earthtime.sav'

 
  noisyspec = poidev(cnoise * time)
  ;noisyspec = cnoise*time
  planet = noisyspec - 2.*(cz + cez + csp + cD + cR + cth)*  time
  sun = (time * cp)/A

  set_colors
  plot, lam, planet/sun, color=0, xrange=[0,1.5]
  oplot, lam, cp/sun*time, color=90


  snr = cp*time/sqrt((cp+2*cb)*time)
  noisy = randomn(seed, N_elements(A))*A/SNR+A

  plot, lam, planet, color=0, xrange=[0,1.5]
  oplot, lam, cp*time, color=90

; wholefeature = where(lam gt 0.84 and lam lt 1.08)
; continuum1 = where(lam gt 0.84 and lam lt 0.88)
; continuum2 = where(lam gt 1 and lam lt 1.08)
; water = where(lam ge 0.88 and lam le 1)
;  plot, lam[wholefeature], planet[wholefeature], color=0, xrange=[0,1.5]
;  oplot, lam[wholefeature], cp[wholefeature]*time, color=90
;  oplot, lam[continuum1], planet[continuum1], color=230
;  oplot, lam[continuum2], planet[continuum2], color=230
  

;  continuum = [lam[continuum1], lam[continuum2]]
;  flx = [planet[continuum1], planet[continuum2]]
;  test1 = spline(continuum, flx, lam[wholefeature],0.1)
;  test2 = poly_fit(continuum, flx, 2)
;  pol = poly(lam[wholefeature], test2)
;oplot, lam[wholefeature], pol, color=100

;whole = planet[wholefeature]/pol
;plot, lam[wholefeature], whole, color=0
;waterline = planet[water]/pol

;equivalent width
;equ = tsum(lam[wholefeature], ((1.-whole)/1.))

;snr of detection
;difference of noiseless spectrum and difference of noise 
;signal = equivalent width
;noise = equivalent width of integrated noise 
;error propegation for equivalent width calculation?
;significance of detection is not the same as signal to noise
;push error array through the equation in email from Aki
;that will give significance of detection of the line

;ummm...is this right??  I have doubts....
;allnoise = (cb[wholefeature]*time)/pol
;equnoise = tsum(lam[wholefeature], ((1.-allnoise)/1.))
 
;print, equ


spawn, 'mkdir plots/'+ whichplanet

pfile = 'plots/'

set_plot, 'ps'
device, file=pfile+whichplanet+'/spec_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, planet, color=0, thick=5, xtitle='wavelength [microns]', ytitle='photon counts'
oplot, lam, cp*time, color=90, thick=5
device, /close
set_plot, 'x'

set_plot, 'ps'
device, file=pfile+whichplanet+'/nonoise_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, cp*time, color=0, thick=5, xtitle='wavelength [microns]', ytitle='photon counts'
device, /close
set_plot, 'x'

set_plot, 'ps'
device, file=pfile+whichplanet+'/noisy_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, planet, color=0, thick=5, xtitle='wavelength [microns]', ytitle='photon counts'
device, /close
set_plot, 'x'

set_plot, 'ps'
device, file=pfile+whichplanet+'/inttime_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, dtsnr, color=0, thick=5, xtitle='wavelength [microns]', ytitle='integration time [hours]', /ylog
device, /close
set_plot, 'x'

;---albedos
set_plot, 'ps'
device, file=pfile+whichplanet+'/albspec_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, planet/sun, color=0, thick=5, xtitle='wavelength [microns]', ytitle='albedo'
oplot, lam, cp*time/sun, color=90, thick=5
device, /close
set_plot, 'x'

set_plot, 'ps'
device, file=pfile+whichplanet+'/albnonoise_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, cp*time/sun, color=0, thick=5, xtitle='wavelength [microns]', ytitle='albedo'
device, /close
set_plot, 'x'

set_plot, 'ps'
device, file=pfile+whichplanet+'/albnoisy_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, planet/sun, color=0, thick=5, xtitle='wavelength [microns]', ytitle='albedo'
device, /close
set_plot, 'x'


;----------noises-----------


set_plot, 'ps'
device, file=pfile+whichplanet+'/noises_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.eps', /color, xsize='8', ysize='5', /encapsulated, bits=8, /inches
plot, lam, cp, color=0, thick=5, xtitle='wavelength [microns]', ytitle='counts', yrange=[0.000001, 100], /ylog
oplot, lam, cz, color=20
oplot, lam, cez, color=40
oplot, lam, csp, color=70
oplot, lam, cd, color=100
oplot, lam, cr, color=120
oplot, lam, cth, color=180

legend, ['planet counts', 'zodi', 'exozodi', 'speckle', 'dark', 'read', 'thermal'], color=[0,20,40,70,100,120,180], linestyle=[0,0,0,0,0,0,0], textcolor=[0,0,0,0,0,0,0,0]

device, /close
set_plot, 'x'


save, lam, dtsnr, planet, cp, time, sun, filename=pfile+whichplanet+'/vars_'+whichplanet+'_'+strtrim(diam,1)+'diam_'+strtrim(res,1)+'res_'+strtrim(wantsnr,1)+'_'+strtrim(lammin,1)+'_'+strtrim(lammax,1)+'_'+strtrim(tput,1)+'tput_'+strtrim(Tsys,1)+'Temp.sav'
stop

  ;stop
END
