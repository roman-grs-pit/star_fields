import numpy as np

def galext_r(X,extinf,b,hred,Zsun):
    '''
                                ; This is the NCI (New COBE/IRAS; Chen
                                ; et al. 1999) extinction model. It
                                ; appears to really break down at b <
                                ; 2.5 deg.  
  
    ; extinf is extinction at infinity, probably derived from galactic_ext.pro
    ; b is latitude in radians (Pi/2 though -Pi/2).
    ; hred is scale height of absorbing material, typically around 100 pc
    ; Zsun is Zsun is our sun's offset from the plane. Assumed positive. Chen et al. (1999) gives a value of 27.5 pc.
    '''

    const1=(1-2*np.exp(Zsun/hred))
    asinb=abs(np.sin(b))
    
    if b >= 0: 
        YMOD=extinf*(1-np.exp(-X*np.sin(b)/hred))
    else: 
        sel = X*asinb <= Zsun
        YMOD = np.zeros(len(X))
        YMOD[sel] = extinf*((1-np.exp(-X[sel]*np.sin(b)/hred))/const1)
        YMOD[~sel] = ((const1+np.exp((2*Zsun+X[~sel]*np.sin(b))/hred))/const1)
    #;   testdoug=[(((1-exp(-X*sin(lat)/hred))/const1)*(X*asinb le Zsun)),((const1+exp((2*Zsun+X*sin(b))/hred))/const1)*(X*asinb gt Zsun)]
    #endelse
    return YMOD
