PRO density_r,  X, P, YMOD

     ; Based on model of Juric et al. (2008), mostly using the
     ; parameters of Chang et al. (2010)
     ; Does not include any Bulge contribution
   
	R0=P(0)     ; Distance from sun to galaxy center (pc). 8000 pc typical (Reid 1993)		 ; 
	Z0=P(1)     ; Distance of the sun from the midplane (pc) 27.5 from Chen et al. (1999). 25 from Juric et al. (2008). Difficult to distiguish from 0.
     ;   n0=P(2)     ; Local thin disk number density of stars. Set to 1 for density normalized to local.
        b=P(2)      ; Galactic Latitude in radians (Pi/2 through -Pi/2).
	long=P(3)  ; Galactic Longitude in radians (0 through 2Pi, 0 pointed at Gal. center).
        n0=1.0
        
	aZ0=abs(Z0)
	Z=X*sin(b)+Z0
	R=sqrt((X*cos(b))^2-2*R0*X*cos(b)*cos(long)+R0^2)

    
        hz1=300.0	  	; Thin disk scale height
        hr1=2600.0		; Thin disk scale length
        f1=1.0			; Thin disk normalization 
        hz2=900.0		; Thick disk scale height
        hr2=3600.0		; Thick disk scale length
	f2=0.12			; Thick disk normalization
        p=2.8			; Power index (spheroidal halo)
        c=0.64			; Flattening parameter (spheroidal halo)
	fh=0.005		; Spheroid normalization


	D1=f1*exp(-((R-R0)/hr1) - ((abs(Z)-aZ0)/hz1))
	D2=f2*exp(-((R-R0)/hr2) - ((abs(Z)-aZ0)/hz2))
	Hs=fh*((R^2+(Z/c)^2)/(R0^2+ (Z0/c)^2))^(-p/2)

	; Return the sum total of density and each component separately, as well.
        YMOD=n0*(D1+D2+Hs)/(1+f2+fh)

END
	
