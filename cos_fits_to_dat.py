import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.io.fits as fits
from astropy.convolution import convolve, Box1DKernel
import glob

"""
Extracts COS spectra from x1dsum.fits files and saves it as an ascii .dat file with columns wavelength, flux, error, dq

Usage:
    
call the function cos_fits_to_dat.fits_to_dat(files=[], plot=True, file_path=os.getcwd()+'/', save_path=os.getcwd()+'/spectra/', filename='long')

If no arguments are called, it will extract all x1dsum.fits files in the working directory, make a new directory 'spectra' in the working directory and save the .dat files there.
    
Arguments:

- files: Array, list of x1dsum.fits files to stitch. Default = []

- plot: boolean. Make a plot of the stitched spectrum at the end. Default = True

- file_path: string, path to where the xld/xlf files are. Remember to include a '/' at the end. Default= working directory

- save_path: string, directory to save the stitched spectra in. Remember to include a '/' at the end. Default = new directory 'stitched_spectra' in the working directory

- file_name: string, what name to save the .dat file to. Two options built in:
    - 'long' (default) : uses header KEYWORDS to save the file as 'TARGNAME_'INSTRUME_DETECTOR_OPT_ELEM_DATE-OBS:TIME-OBS_ROOTNAME.dat'
    - 'short' : saves the file as 'ROOTNAME.dat'

"""

def filewriter(wavelength, flux, error, dq, save_path, filename):
    """
    write the spectrum to file in savepath/file name in space-separated columns of wavelength, flux, flux_error, dq.
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)    
    fl=open((save_path+filename),'w') 
    for w, f, e, q in zip(wavelength, flux, error, dq):
        fl.write('%f %g %g %i\n'%(w,f,e,q))
        
def plot_spectrum(wavelength, flux, error, dq, rootname):
    """diagnostic plot. Smooths with a 5pt boxcar. Plots the raw spectrum in grey, then overlays another spectrum with flagged points removed"""
	flux = convolve(flux,Box1DKernel(5))
	plt.figure(rootname)
	plt.plot(wavelength, flux, '0.5')
	wavelength_dq, flux_dq = wavelength[dq==0], flux[dq==0]
	plt.plot(wavelength_dq, flux_dq, label = 'spectrum')
	plt.plot(wavelength, error, label='error')  
	plt.legend()
	plt.show()

def get_data(data):
    """extract the data from the fits file and combine them into numpy arrays"""
    w = np.array([], dtype=float)
    f = np.array([], dtype=float)
    e = np.array([], dtype=float)
    dq = np.array([], dtype=int)
    
    for j in range(len(data)):
        w = np.concatenate((w, data['wavelength'][j]))
        f = np.concatenate((f, data['flux'][j]))
        e = np.concatenate((e, data['error'][j]))
        dq = np.concatenate((dq, data['dq'][j]))

    arr1inds = w.argsort()# reorders the various segements by wavelength
    w = w[arr1inds]
    f = f[arr1inds]
    e = e[arr1inds]
    dq = dq[arr1inds]

    return w, f, e, dq 

def fits_to_dat(files=[], plot=True, file_path=os.getcwd()+'/', save_path=os.getcwd()+'/spectra/', filename='long'):
    #main function, gathers fits files
    if files == []:
        fitsfiles = glob.glob(file_path+'*x1dsum.*')  
   

    for fitsfile in fitsfiles:
        hdul=fits.open(fitsfile)
        hdr = hdul[0].header
        hdr2 = hdul[1].header
        data = hdul[1].data
        hdul.close()

        w, f, e, dq = get_data(data)

        #write to file
        if filename == 'long':
            filename = hdr['TARGNAME']+'_'+hdr['INSTRUME']+'_'+hdr['DETECTOR']+'_'+hdr['OPT_ELEM']+'_'+hdr2['DATE-OBS']+':'+hdr2['TIME-OBS']+'_'+hdr['ROOTNAME']+'.dat'
        if filename == 'short':
            filename = hdr['ROOTNAME']+'.dat'
        filewriter(w, f, e, dq, save_path = save_path, filename = filename  )

        #make a plot
        if plot == True:
            plot_spectrum(w, f, e, dq, filename)

            