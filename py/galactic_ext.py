import healpy as hp
import os
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord


def galactic_ext(RA,DEC,resol='12.0', wavelength=1):
    '''
    RA and Dec expected in degrees in FK5 J2000 system

    Wavelength is desired wavelength of output extinction, A_lambda, and must be input in microns. 
        A_lambda/A_J assumed proportional to a power law with a power given as:  power=1.7     
        Bailey, M.E., and Williams, D.A., eds. 1988. Dust in the Universe, Cambridge: Cambridge University Press, 573
            Should work well for Y,J,H,K and be
            ok for L,M,R,I 
    resol: Pick resolution of gaussian beam used to derive extinction: 1) 3.0" 2) 4.5" 3) 12.0"
    outputs expected Galactic extinction
    '''
    allskyext_fn=os.getenv('allskyext')+'/NICER_AJ_M2a_FWHM'+resol+'.fits'#'/Users/colbert/AllSkyExtinction/'
    if not os.path.isfile(allskyext_fn):
        print('ERROR extinction file '+allskyext_fn +'not found; maybe $allskyext environment variable is not set correctly?')
                 
    extmap = hp.read_map(allskyext_fn)
    rd = SkyCoord(ra=RA*u.degree,dec=DEC*u.degree, unit='deg', frame='fk5')
    gc = rd.transform_to('galactic')
    long = gc.l.deg
    lat = gc.b.deg
    phi = np.pi/2*(1-lat/90.)
    th = 2*np.pi*(long/360)
    pix = hp.ang2pix(2048,phi,th,nest=True)
    extval = extmap[pix]
    power = 1.7
    output = extval*(1.25/wavelength)**power
    return output
