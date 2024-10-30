pro starcount, RA, DEC, wavelength, type, outmags, outdensity

    allskyext='/Users/colbert/AllSkyExtinction/'
    numpc=5e4  ; Number of steps to integrate over
    dr=1.0	   ; Each step in terms of pc
    R_gal=8000.0  ; Distance from sun to galaxy center (pc)
    scaleh=100.0   ;Scale height of dust disk (pc)
    ZSol=25.0	   ; Solar offset from galactic plane (pc)

                                ; Type indicates whether input is expected in RA/Dec (deg): type=1
                                ; or Gal. longitude and latitude: type=2 
    
    ; RA and Dec expected in degrees in FK5 J2000 system
    if type eq 1 then begin
      glactc, RA,DEC,2000.0,long1,lat1,1,/degree
      print,'Longitude: ',long1,'  Latitude:',lat1
    endif
    if type eq 2 then begin
      glactc, RA1,DEC1,2000.0,RA,Dec,2,/degree
      lat1=RA
      long1=DEC
      print,'Right Ascension: ',RA1,'  Declination:',DEC1
      RA=RA1
      DEC=DEC1
    endif
     
     ; Wavelength is desired wavelength of output extinction, A_lambda, and
     ;    must be input in microns.
     galactic_ext,3,RA,DEC,wavelength,galext0
     if galext0 lt 0.0 then begin
         print,"Extinction less than 0, set to zero."
         galext0=0.0
     endif

      ; Convert Galactic coords to radians
      lat=(!pi/2.0)*(lat1/90.0)
      long=(2*!pi)*(long1/360.0)
      print,'Galactic Extinction at Infinity: ',galext0
      
      ; Derive galactic extinction as function of distance from Sun
      Radius=(findgen(numpc-1)+1)*dr+1
      params1=[galext0,lat,scaleh,ZSol]
      galext_r,Radius,params1,galext
     
      ; Derive density of stars (relative to local) as function of distance from the Sun
      params2=[R_gal,ZSol,lat,long]
      density_r,Radius,params2,density
 
      ;Read in stellar luminosity function  
      readcol,allskyext+'stellarlumfunction.dat',absK,rho1,rho1E,rho2,format='f,f,f,f'
      readcol,allskyext+'stellarHK.dat',absK,HKcolor,HKcolorGIANT
      modifier=1/33510.0   ; Account for fact density of lum funct. given per 20pc radius sphere

      deltamag=0.5
      numbins=(31+10)/deltamag
      num_extinct=3
      lth=n_elements(absK)
      density_all=dblarr(lth,num_extinct,numbins,2)
      allmags=float(findgen(numbins))*deltamag-10
      
      for j=0,numpc-2 do begin
         R=Radius(j)
         for i=0,lth-1 do begin           
                 m=absK(i)+5*alog10(R)-5+galext(j)+HKcolor(i)+1.39 ; Last component switches to AB. 
                 count_slice=4*!pi*modifier*rho2(i)*density(j)*(R^2)*dr
                 m2=absK(i)+5*alog10(R)-5+galext(j)+HKcolorGIANT(i)+1.39
                 count_slice2=4*!pi*modifier*(rho1(i)-rho2(i))*density(j)*(R^2)*dr
                 for k=0,numbins-1 do begin
                    if m gt allmags(k)-deltamag/2.0 and m lt allmags(k)+deltamag/2.0 then begin  
                       if galext(j) le 0.25 then density_all(i,0,k,0)=density_all(i,0,k,0)+count_slice*(1/41253.0) ; Density of this Abs Mag per square degree, M.S. only
                       if (galext(j) gt 0.25) AND (galext(j) le 0.75) then density_all(i,1,k,0)=density_all(i,1,k,0)+count_slice*(1/41253.0)
                       if (galext(j) gt 0.75) AND (galext(j) le 1.25) then density_all(i,2,k,0)=density_all(i,2,k,0)+count_slice*(1/41253.0)
                       if galext(j) gt 1.25 then print, 'Warning: Extinction greather than 1.25 mags. Not included in counts. MS. '
                    endif
                    if  m2 gt allmags(k)-deltamag/2.0 and m2 lt allmags(k)+deltamag/2.0 then begin
                       if galext(j) le 0.25 then density_all(i,0,k,1)=density_all(i,0,k,1)+count_slice2*(1/41253.0) ; Density of this Abs Mag per square degree,only Giants or W.D.
                       if (galext(j) gt 0.25) AND (galext(j) le 0.75) then density_all(i,1,k,1)=density_all(i,1,k,1)+count_slice2*(1/41253.0)
                       if (galext(j) gt 0.75) AND (galext(j) le 1.25) then density_all(i,2,k,1)=density_all(i,2,k,1)+count_slice2*(1/41253.0)
                       if galext(j) gt 1.25 then print, 'Warning: Extinction greather than 1.25 mags. Not included in counts. GIANTS. '
                    endif
                 endfor
                 
	;	 print,absK(i),R,m,count_slice*(1/41253.0)
	  endfor
      endfor

      density_combined=fltarr(numbins,2)
      for k=0,numbins-1 do begin
         for i=0,lth-1 do begin
             density_combined(k,0)=density_combined(k,0)+density_all(i,0,k,0)
             density_combined(k,1)=density_combined(k,1)+density_all(i,0,k,1)
          endfor
      endfor
      outdensity=density_all
      outmags=allmags
  ;   forprint,allmags,density_combined(*,0),density_all(25,*,0),density_all(20,*,0),density_all(15,*,0),density_all(10,*,0),density_all(5,*,0),density_all(0,*,0)
  ;   bob=where(allmags lt 10.0 and allmags gt 8)
  ;   print,total(density_combined(bob,0))*0.0156
     
end    
