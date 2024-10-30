PRO galext_r, X, P, YMOD
                                ; This is the NCI (New COBE/IRAS; Chen
                                ; et al. 1999) extinction model. It
                                ; appears to really break down at b <
                                ; 2.5 deg.  
  
    ; P(0) is extinction at infinity, probably derived from galactic_ext.pro
    ; P(1) is latitude in radians (Pi/2 though -Pi/2).
    ; P(2) is scale height of absorbing material, typically around 100 pc
    ; P(3) is Zsun is our sun's offset from the plane. Assumed positive. Chen et al. (1999) gives a value of 27.5 pc.

    galext0=P(0)
    b=P(1)
    hred=P(2)
    Zsun=P(3)

    const1=(1-2*exp(Zsun/hred))
    asinb=abs(sin(b))
    
    if b ge 0 then begin
       YMOD=galext0*(1-exp(-X*sin(b)/hred))
    endif else begin
       YMOD=galext0*(((1-exp(-X*sin(b)/hred))/const1)*(X*asinb le Zsun)+((const1+exp((2*Zsun+X*sin(b))/hred))/const1)*(X*asinb gt Zsun))
    ;   testdoug=[(((1-exp(-X*sin(b)/hred))/const1)*(X*asinb le Zsun)),((const1+exp((2*Zsun+X*sin(b))/hred))/const1)*(X*asinb gt Zsun)]
    endelse
END
