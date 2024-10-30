#pro starcount, RA, DEC, wavelength, type, outmags, outdensity

import os
import argparse

import astropy.units as u
from astropy.coordinates import SkyCoord

from density_r import density_r
from galext_r import galext_r
from galactic_ext import galactic_ext

import numpy as np

#parser = argparse.ArgumentParser()
#parser.add_argument('--type',help=' Type indicates whether input is expected in RA/Dec (deg): type=1,or Gal. longitude and latitude: type=2 ', options=[1,2],default=1,dtype=int)
#parser.add_argument('--RA',help='RA and Dec expected in degrees in FK5 J2000 system', required=True)
#parser.add_argument('--DEC',help='RA and Dec expected in degrees in FK5 J2000 system', required=True)
#parser.add_argument('--wavelength',help='Wavelength is desired wavelength of output extinction, A_lambda, and must be input in microns', default=1)
#args = parser.parse_args()
numpc=5e4  # Number of steps to integrate over
dr=1.0     # Each step in terms of pc
R_gal=8000.0  # Distance from sun to galaxy center (pc)
scaleh=100.0   #Scale height of dust disk (pc)
ZSol=25.0      # Solar offset from galactic plane (pc)


allskyext=os.getenv('allskyext')
star_data=os.getenv('star_data')

def get_starcounts(tp,RA,DEC,wavelength):


    if tp == 1:
        #RA and Dec expected in degrees in FK5 J2000 system
        rd = SkyCoord(ra=RA*u.degree,dec=DEC*u.degree, unit='deg', frame='fk5')
        gc = rd.transform_to('galactic')
        long1 = gc.l.deg
        lat1 = gc.b.deg
        #RA = args.RA
        #DEC = args.DEC
        print('Longitude: '+str(long1),'Latitude '+str(lat1))
    elif tp == 2:
        gc = SkyCoord(l=RA*u.degree, b=DEC*u.degree, frame='galactic')
        rd = gc.transform_to('fk5')
        lat1 = args.RA
        long1 = args.DEC
        RA = rd.ra.deg
        DEC = rd.dec.deg
        print('Right Ascension: '+str(RA),'Declination: '+str(DEC))
    else:
        sys.exit('Unsupported option for '+str(args.type))
    
    galext0 = galactic_ext(RA,DEC,resol='12.0', wavelength=wavelength)
    if galext0 < 0:
       galext0 = 0
         
    #Convert Galactic coords to radians
    lat=(np.pi/2.0)*(lat1/90.0)
    long=(2*np.pi)*(long1/360.0)
    print('Galactic Extinction at Infinity: ',galext0)
          
    #Derive galactic extinction as function of distance from Sun
    
    #Radius=(findgen(numpc-1)+1)*dr+1
    Radius = (np.arange(numpc-1)+1)*dr+1
    #params1=[galext0,lat,scaleh,ZSol]
    gext_r = galext_r(Radius,galext0,lat,scaleh,ZSol)
    #galext_r,Radius,params1,galext
         
    #Derive density of stars (relative to local) as function of distance from the Sun
    params2=[R_gal,ZSol,lat,long]
    dens_r = density_r(Radius,R_gal,ZSol,lat,long)
    #density_r,Radius,params2,density
     
    #Read in stellar luminosity function  
    lum_fn = star_data+'/stellarlumfunction.dat'
    lum_dat = np.loadtxt(lum_fn).transpose()
    absK = lum_dat[0]
    rho1 = lum_dat[1]
    rho1E = lum_dat[2]
    rho2 = lum_dat[3]
          
    #      readcol,allskyext+'stellarlumfunction.dat',absK,rho1,rho1E,rho2,format='f,f,f,f'
    HK_fn = star_data+'/stellarHK.dat'
    HK_dat = np.loadtxt(HK_fn).transpose()
    absK = HK_dat[0]
    HKcolor = HK_dat[1]
    HKcolorGIANT = HK_dat[2]
    #      readcol,allskyext+'stellarHK.dat',absK,HKcolor,HKcolorGIANT
    modifier=1/33510.0   # Account for fact density of lum funct. given per 20pc radius sphere
    
    deltamag=0.5
    numbins=int((31+10)/deltamag)
    num_extinct=3
    lth=len(absK)
    density_all= np.zeros((lth,num_extinct,numbins,2))
    #      density_all=dblarr(lth,num_extinct,numbins,2)
    allmags = np.arange(numbins)*deltamag-10
    #      allmags=float(findgen(numbins))*deltamag-10
          
    for j in range(0,int(numpc)-1): 
        R=Radius[j]
        for i in range(0,lth):            
            m=absK[i]+5*np.log10(R)-5+gext_r[j]+HKcolor[i]+1.39 #; Last component switches to AB. 
            count_slice=4*np.pi*modifier*rho2[i]*dens_r[j]*(R**2)*dr
            m2=absK[i]+5*np.log10(R)-5+gext_r[j]+HKcolorGIANT[i]+1.39
            count_slice2=4*np.pi*modifier*(rho1[i]-rho2[i])*dens_r[j]*(R**2)*dr
            for k in range(0,numbins): 
                if m > allmags[k]-deltamag/2.0 and m < allmags[k]+deltamag/2.0:
                    if gext_r[j] <= 0.25:
                        density_all[i][0][k][0]=density_all[i][0][k][0]+count_slice*(1/41253.0) #; Density of this Abs Mag per square degree, M.S. only
                    if (gext_r[j] > 0.25) and (gext_r[j] <= 0.75):
                        density_all[i][1][k][0] =density_all[i][1][k][0]+count_slice*(1/41253.0)
                    if (gext_r[j] > 0.75) and (gext_r[j] <= 1.25): 
                        density_all[i][2][k][0]=density_all[i][2][k][0]+count_slice*(1/41253.0)
                    if gext_r[j] > 1.25:
                        print('Warning: Extinction greather than 1.25 mags. Not included in counts. MS. ')
                       
                if  m2 > allmags[k]-deltamag/2.0 and m2 < allmags[k]+deltamag/2.0: 
                    if gext_r[j] <= 0.25:
                        density_all[i][0][k][1]=density_all[i][0][k][1]+count_slice2*(1/41253.0) #; Density of this Abs Mag per square degree,only Giants or W.D.
                    
                    if (gext_r[j] > 0.25) and (gext_r[j] <= 0.75):
                        density_all[i][1][k][1]=density_all[i][1][k][1]+count_slice2*(1/41253.0)
                    if (gext_r[j] > 0.75) and (gext_r[j] <= 1.25):
                        density_all[i][2][k][1] = density_all[i][2][k][1]+count_slice2*(1/41253.0)
                    if gext_r[j] > 1.25:
                        print('Warning: Extinction greather than 1.25 mags. Not included in counts. GIANTS. ')
    
                     
            #print(absK[i],R,m,count_slice*(1/41253.0))
    
    density_combined=np.zeros((numbins,2))
    for k in range(0,numbins): 
        for i in range(0,lth): 
            density_combined[k][0]=density_combined[k][0]+density_all[i][0][k][0]
            density_combined[k][1]=density_combined[k][1]+density_all[i][0][k][1]
    outdensity=density_all
    outmags=allmags
    return outdensity,outmags
    #  ;   forprint,allmags,density_combined(*,0),density_all(25,*,0),density_all(20,*,0),density_all(15,*,0),density_all(10,*,0),density_all(5,*,0),density_all(0,*,0)
    #  ;   bob=where(allmags lt 10.0 and allmags gt 8)
    #  ;   print,total(density_combined(bob,0))*0.0156
         
    #end    
