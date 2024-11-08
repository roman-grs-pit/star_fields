import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
# Spectra tools
import pysynphot as S
import webbpsf

import grizli
import grizli.utils

roman_base_dir = os.getenv('HOME')+'/Dropbox/RomanGRS/products/' #relying on 2022 sim products to start with
direct_fits = roman_base_dir+"FOV0/roll_0/dither_0x_0y/SCA1/GRS_FOV0_roll0_dx0_dy0_SCA1_direct_final.fits"
ready_fits = roman_base_dir+"FOV0/roll_0/dither_0x_0y/SCA1/ready_direct_GRS_FOV0_roll0_dx0_dy0_SCA1_direct_final.fits"
empty_seg = roman_base_dir+"FOV0/roll_0/dither_0x_0y/SCA1/empty_seg.fits"

file = fits.open(direct_fits)

file[0].header["INSTRUME"] = "ROMAN"
file[0].header["FILTER"] = "d1_"
file[1].header["CONFFILE"] = roman_base_dir+"configuration/Roman.det1.07242020.conf" # This had to be a path, not just a filename; otherwise, grizli can't find the sensitivity fits

file.writeto(ready_fits, overwrite=True)

header = file[1].header # We need to copy this header info into the segmentation file

file.close()

SED_dir = roman_base_dir+"FOV0/SEDs/" # Change to your path to directory containing SEDs

# Create F158 Filter Bandpass object
df = Table.read(os.path.join(SED_dir, "wfirst_wfi_f158_001_syn.fits"), format='fits')
bp = S.ArrayBandpass(df["WAVELENGTH"], df["THROUGHPUT"])

from grizli.model import GrismFLT

pad = 100 # Objects near the edge of the detector may have a trace with a center of the detector. Padding ensures it will still disperse

#roman = GrismFLT(direct_file=ready_fits, seg_file=empty_seg, pad=pad)
roman_temp = GrismFLT(direct_file=ready_fits, seg_file=empty_seg, pad=pad)
roman_flat = GrismFLT(direct_file=ready_fits, seg_file=empty_seg, pad=pad)

# STAR generic circle shape
gen_circle = lambda x, y, x_0, y_0: (x-x_0)**2 + (y-y_0)**2

def get1d_spec(image,minl=1e4,maxl=2e4,res=25,minx=-400,maxx=599):
    #extract 1D spectrum from grism image
    #sums trace based on dispersion model applied to segment of direct image
    #~12 angstroms per pixel
    nbin = int((maxl-minl)*1.0001/res)
    fluxl = np.zeros(nbin)
    npixl = np.zeros(nbin)
    segv = image.seg.nonzero()
    lamvl = np.arange(minl+res/2.,maxl,res)
    for i in range(0,len(segv[0])):
        xp = segv[1][i]
        yp = segv[0][i]
        for dx1 in range(minx,maxx):
            dy, lam = image.conf.get_beam_trace(x=xp-image.pad[0], y=yp-image.pad[1], dx=dx1, beam='A')
            dyv = dy
            dxv = dx1
            lamv = lam
            if lamv > minl and lamv < maxl:
                ypf = int(yp+dyv)
                xpf = int(xp+dxv)
                if ypf < len(image.model[0]) and xpf < len(image.model[0]):
                    modv = image.model[ypf][xpf]
                    binv = int((lamv-minl)/res)
                    fluxl[binv] += modv
                    npixl[binv] += 1.
    return lamvl,fluxl,npixl

def plot_extract_star_spec(star_type='G0V',mag=22,cat="phoenix",github_dir='/Users/ross.1333/Documents/GitHub/',minlam=1e4,maxlam=2e4):
    N = 0
    file = fits.open(direct_fits)
    # Read in star spectrum
    if cat == "phoenix":
        src = webbpsf.specFromSpectralType(star_type, catalog="phoenix")
        wave = np.linspace(10000, 20000, 10000)
        flux = src(wave).value
    if cat == "JC":
        tempdir = github_dir+'star_fields/data/SEDtemplates/'
        temp = np.loadtxt(tempdir+star_type).transpose()
        wave = temp[0]
        sel = wave > minlam
        sel &= wave < maxlam
        wave = wave[sel]
        flux = temp[1]
        flux = flux[sel]
    star_spec = S.ArraySpectrum(wave=wave, flux=flux, waveunits="angstroms", fluxunits="flam")
    
    temp_seg = np.zeros((4288,4288), dtype="float32") # Grizli expects arrays to have dtype="float32"
    idu= 1
    x_0 = 2000 #row["X_IMAGE"] + 100 # Add 100 to account for padding
    y_0 = 2000 #row["Y_IMAGE"] + 100 
    
    radius = 5 
    circ = lambda x, y: gen_circle(x, y, x_0, y_0) <= (radius ** 2) # faster to compare squares than sqrt (I think)
        
    # Define cutout bounds
    x_min = max(int(x_0 - radius + 1), 0)       # max/min ensures we stay on the detector (plus padding)
    x_max = min(int(x_0 + radius + 1), 4288)    # The +1 is so we've got a bit of padding, to ensure our cutout's not rounding down too far
    y_min = max(int(y_0 - radius + 1), 0)
    y_max = min(int(y_0 + radius + 1), 4288)
        
    # Create Meshgrid within bounds
    x = np.arange(x_min, x_max)
    y = np.arange(y_min, y_max)
    x_grid, y_grid = np.meshgrid(x, y)
    
    condition = circ(x_grid, y_grid)
    
    temp_seg[y_min:y_max, x_min:x_max][np.where(condition)] = idu #

    roman_temp.seg = temp_seg
    roman_flat.seg = temp_seg
    # SIMULATION

    # Adjust spectrum scale using F158 Filter Bandpass object and star magnitude
    spec = star_spec.renorm(mag, "abmag", bp)
    spec.convert("flam")

    # By default, grizli tries to compute a cutout size. This cutout size is not large enough for the roman grism.
    # In 4) FOV0_sims/notebooks/dy-by-optimize.ipynb, Keith estimates the maximum needed size to be 77 for detector 1.
    # See that notebook for more details
    roman_temp.compute_model_orders(id=idu, mag=mag, compute_size=False, size=77, in_place=True, store=False,
                       is_cgs=True, spectrum_1d=[spec.wave, spec.flux])

    lamv,flux,pix = get1d_spec(roman_temp)

    flatflux = np.ones(len(wave))
    flatspec = S.ArraySpectrum(wave=wave, flux=flatflux, waveunits="angstrom", fluxunits="flam").redshift(0)
    flatspec = flatspec.renorm(mag, "abmag", bp)
    flatspec.convert("flam")
    roman_flat.compute_model_orders(id=idu, mag=mag, compute_size=False, size=77, in_place=True, store=False,
                           is_cgs=True, spectrum_1d=[flatspec.wave, np.ones(len(flatspec.wave))])#flatspec.flux])
    _,flatflux,flatpix = get1d_spec(roman_flat)

    plt.plot(lamv,flux/flatflux,label='simulated extracted spectrum')
    plt.plot(spec.wave, spec.flux,label='input',zorder=0)

    plt.title('stellar type '+star_type+' F1500 mag '+str(round(mag,1)))
    plt.ylabel('flux')
    plt.xlabel('wavelength (angstroms)')
    plt.legend()
    #plt.show()        

if __name__ == '__main__':
    github_dir='/Users/ross.1333/Documents/GitHub/'
    figs = []
    #types in pheonix library
    stellar_types = ['F0V','K0V','M0V','B0V','G0V']
    for tp in stellar_types:
        fig = plt.figure()
        plot_extract_star_spec(star_type=tp)
        figs.append(fig)
    templates = open(github_dir+'star_fields/data/SEDtemplates/input_spectral_STARS.lis').readlines()
    tmpl = []
    for i in range(0,len(templates)):#ind in temp_ind:
        tmp = templates[i].strip('\n')
        if 'garbage' not in tmp:
            tmpl.append(tmp)
        #if tmp not in tmpl:
        #    tmpl.append(tmp.strip('\n'))


    for tmp in tmpl:
        fig = plt.figure()
        plot_extract_star_spec(star_type=tmp,cat='JC')
        figs.append(fig)

    with PdfPages('plots/cannonical_star_spec.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
