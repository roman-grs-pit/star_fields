pro galactic_ext, resol, RA, DEC, wavelength, output

		  allskyext='/Users/colbert/AllSkyExtinction/'
                  ; RA and Dec expected in degrees in FK5 J2000 system
  
                  ; Wavelength is desired wavelength of output extinction, A_lambda, and
                  ; must be input in microns.
                  ; A_lambda/A_J assumed proportional to a power law with a
                  ; power given as: 

                  power=1.7     ; Bailey, M.E., and Williams, D.A., eds. 1988. Dust in the Universe, Cambridge: Cambridge University Press, 573
                                ; Should work well for Y,J,H,K and be
                                ; ok for L,M,R,I 
                  
    		  ; Pick resolution of gaussian beam used to derive extinction: 1) 3.0" 2) 4.5" 3) 12.0"
    		  case resol of
		       1:  input=allskyext+'NICER_AJ_M2a_FWHM3.0.fits'
		       2:  input=allskyext+'NICER_AJ_M2a_FWHM4.5.fits'
		       3:  input=allskyext+'NICER_AJ_M2a_FWHM12.0.fits'
		  endcase
                 
		 ; Convert RA+Dec to Galactic coords, then covert to radians
		  glactc, ra,dec,2000.0,long,lat,1,/degree
		  lat2=(!pi/2.0)*(1-lat/90.0)
		  long2=(2*!pi)*(long/360.0)
                ;  forprint,long,lat,long2,lat2
                  
                 ; Using HEALPix routine, convert
                                ; lat+long to nested pixel coordinates.
                  ang2pix_nest,2048.0,lat2,long2,values

                  ; Read in NICER map (Juvela & Montillaud 2016) of J-Band
                  ; extinction, using HEALPix IDL routine, read_fits_map
		  read_fits_map,input,data,head,ehead
                  temp=data(values)
                  output=temp*(1.25/wavelength)^(power)
end
		  
		  
