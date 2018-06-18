#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate 
import astropy.io.fits as fits
import glob

"""
Stitches STIS echelle spectra with overlaping orders. 

At each overlap it interpolates the flux of the order with smaller wavelength grid spacing onto the larger wavelength grid (although the difference is small) and coadds the fluxes.  

Idealy uses x1f files produced using stisblazefix (https://stisblazefix.readthedocs.io). Can use x1d files but these are badly affected be echelle ripple.

For best results with E140M spectra, ask me to do a correction with the updated pht files first. This marginally improves on stisblazefix alone.

Usage:

call stis_echelle_coadd()

If no arguments are called, it will stitch all x1f (or x1d if no x1f files are present) in the working directory, make a new directory in the working directory and save the stiched spectra there. 

Arguments:

- files: Array, list of x1f or x1d files to stitch. Default = []

- plot: boolean. Make a plot of the stiched spectrum at the end. Default = True

- nclip: int, points to clip off each echelle order to clear up the order end problems that are inherent to stis data. 
         I don't know a reason to change it, but the option is there. Default =5 

- file_path: string, path to where the xld/xlf files are. Remember to include a '/' at the end. Default= working directory 

- save_path: string, directory to save the stitched spectra in. Remember to include a '/' at the end. Default = new directory 'stitched_spectra' in the working directory

- file_name: string, what name to save the .dat file to. Two options built in:
    - 'long' (default) : uses header KEYWORDS to save the file as 'TARGNAME_'INSTRUME_DETECTOR_OPT_ELEM_TDATEOBS:TTIMEOBS_ROOTNAME_stitched.dat'
    - 'short' : saves the file as 'ROOTNAME.dat'

-flat: Applies Knox Long's calibration flat to E140M data, currently private.

Output:

A .dat file with space-separated columns of wavelength, flux, flux_error and data quality.

"""

def spectra_adder(fluxes, errors):
    """
    combines the flux at each overlap
    """
    weight_f = np.average(fluxes, axis =0, weights=(1/errors**2))
    weight_e = np.average((weight_f - fluxes)**2, axis=0, weights = (1/errors**2))**0.5
    return weight_f, weight_e

def dq_combine(dq0, dq1):
    """
    makes a new dq array by taking the largest dq value for each point at the overlap. Imperfect but it works, may improve in the future
    """
    dq_comb = []
    for a, b in zip(dq0, dq1):
        dq_comb.append(np.max([a,b]))
    return dq_comb
    
def end_clip(array, nclip):
    """
    clips the first and last nclip pixels off each order
    """
    new_array = []
    shape = np.shape(array)
    order = 0
    while order < shape[0]:
        new_array.append(array[order][nclip:-(nclip+1)])
        order += 1
    return new_array
    
def echelle_coadd(wavelength, flux, err, dq):
    """
    combines echelle orders into one spectrum, stiching them together at the overlap 
    """
    
    #new arrays to put the output in
    w_full = np.array([], dtype=float)
    f_full = np.array([], dtype=float)
    e_full = np.array([], dtype=float)
    dq_full = np.array([], dtype=float)    

    shape = np.shape(flux)
    order = 0
    while order < (shape[0]):
        
        #first add the part that does not overlap ajacent orders to the final spectrum
        if order == 0: #first and last orders do not overlap at both ends
            overmask = (wavelength[order] > wavelength[order + 1][-1])
        elif order == shape[0]-1:
            overmask = (wavelength[order] < wavelength[order - 1][1])
        else:
            overmask = (wavelength[order] > wavelength[order + 1][-1]) & (wavelength[order] < wavelength[order - 1][1])
        w_full = np.concatenate((w_full, wavelength[order][overmask]))
        f_full = np.concatenate((f_full, flux[order][overmask]))
        e_full = np.concatenate((e_full, err[order][overmask]))
        dq_full = np.concatenate((dq_full, dq[order][overmask]))
        
        if order != shape[0]-1:
            
            #interpolate each order onto the one beneath it, with larger wavelength bins. Code adapted from stisblazefix
            f = interpolate.interp1d(wavelength[order + 1], flux[order + 1], fill_value='extrapolate')
            g = interpolate.interp1d(wavelength[order + 1], err[order + 1], fill_value='extrapolate')
            q = interpolate.interp1d(wavelength[order + 1], dq[order + 1], kind='nearest',fill_value=4)
            overlap = np.where(wavelength[order] <= wavelength[order + 1][-1])
            f0 = flux[order][overlap]
            f1 = f(wavelength[order][overlap])
            g0 = err[order][overlap]
            g1 = g(wavelength[order][overlap])
            dq0 = dq[order][overlap]
            dq1 = np.round(q(wavelength[order][overlap])).astype(int)
             
            #combine flux and error at overlap and add to final spectrum
            dq_comb = dq_combine(dq0, dq1)
            w_av = wavelength[order][overlap]
            f_av, e_av = spectra_adder(np.array([f0,f1]),np.array([g0,g1]))
            w_full = np.concatenate((w_full, w_av))
            f_full = np.concatenate((f_full, f_av))
            e_full = np.concatenate((e_full, e_av))
            dq_full = np.concatenate((dq_full, dq_comb))
        order += 1
    
    #stis orders are saved in reverse order, so combined spectra are sorted by the wavelength array
    arr1inds = w_full.argsort()
    sorted_w = w_full[arr1inds]
    sorted_f = f_full[arr1inds]
    sorted_e = e_full[arr1inds]
    sorted_dq = dq_full[arr1inds]
    
    return sorted_w, sorted_f, sorted_e, sorted_dq


def filewriter(wavelength, flux, error, dq, save_path, filename):
    """
    write the spectrum to file in savepath/file name in space-separated columns of wavelength, flux, flux_error, dq.
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)    
    fl=open((save_path+filename),'w') 
    for w, f, e, q in zip(wavelength, flux, error, dq):
        fl.write('%f %g %g %i\n'%(w,f,e,q))
    print('Spectrum saved as '+filename)

def plot_spectrum(wavelength, flux, error, dq, rootname):
    """
    plot the spectrum, error and spectrum with flagged pixels removed
    """
    wavelength, flux, error, dq = wavelength[wavelength > 1160], flux[wavelength > 1160], error[wavelength > 1160], dq[wavelength > 1160] #cut off the scrappy bit at the beginning of E140M spectra
    plt.figure(rootname)
    plt.plot(wavelength, flux, '0.5', label='spectrum')
    plt.plot(wavelength[dq==0], flux[dq==0], label = 'dq filtered spectrum')
    plt.plot(wavelength, error, label='error')  
    plt.legend()
    plt.show()    

def flat_correct(wavelength, flux, error, flatfile):
    """
    applies Knox Long's flat correction to an x1d file. Flat not yet publically available, sorry. 
    """
    print('Applying flat correction')
    np.seterr(divide='ignore', invalid='ignore')
    hdul = fits.open(flatfile)
    flat = hdul[1].data
    hdul.close()
    for i in range(len(wavelength)):
        ft = interpolate.interp1d(flat['WAVELENGTH'][i], flat['FLUX'][i],bounds_error=False, fill_value=1)
        flux[i] /= ft(wavelength[i])
        error[i] /= ft(wavelength[i])
    return flux, error
    
def stis_echelle_coadd(files=[], plot=True, nclip=15, file_path=os.getcwd()+'/', save_path=os.getcwd()+'/stitched_spectra/', flat='', filename='long'):
    """
    main funtion, gather data then stitch, coadd, save and plot them
    """
    
    #get all x1f or x1d files
    if files == []:
        if len(glob.glob(file_path+'*x1d.fits')) > 0 and flat != '':
            print('You are using x1d files with a flat.')
            files = glob.glob(file_path+'*x1d.fits')
        elif len(glob.glob(file_path+'*x1f.fits')) > 0:
            print('You are using x1f files.')
            files = glob.glob(file_path+'*x1f.fits')
        elif len(glob.glob(file_path+'*x1d.fits')) > 0:
            print('You are using x1d files. The result will be better with xlf files made using stisblazefix.')
            files = glob.glob(file_path+'*x1d.fits')
        else:
            print ('There are no x1f or x1d files in file_path :(.')
            os._exit(1)
    

    for fitsfile in files:
        
        #get required data from the x1f/x1d
        print('Working on '+fitsfile)
        hdul = fits.open(fitsfile)
        hdr =  hdul[0].header
        filedata = hdul[1].data
        hdul.close()
        
        #clip the ends off each order
        flux = end_clip(filedata['flux'], nclip)
        dq = end_clip(filedata['dq'], nclip)
        err = end_clip(filedata['error'], nclip)
        wavelength = end_clip(filedata['wavelength'], nclip)
        
        #apply flat correction
        if flat != '' and hdr['OPT_ELEM'] == 'E140M':
            flux, err = flat_correct(wavelength, flux, err, flat)
        
        #stitch the orders together
        w, f, e, dq = echelle_coadd(wavelength, flux, err, dq)
        
        #generate the file name and save the data
        if filename == 'long':
            built_name = hdr['TARGNAME']+'_'+hdr['INSTRUME']+'_'+hdr['DETECTOR']+'_'+hdr['OPT_ELEM']+'_'+hdr['TDATEOBS']+':'+hdr['TTIMEOBS']+'_'+hdr['ROOTNAME']+'_stitched.dat'
        if filename == 'short':
            built_name = hdr['ROOTNAME']+'.dat'
        filewriter(w, f, e, dq, save_path, built_name)
        
        #plot the finished spectrum
        if plot==True:
            plot_spectrum(w, f, e, dq, built_name)


    