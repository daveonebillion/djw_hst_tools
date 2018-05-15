import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.io.fits as fits
from astropy.convolution import convolve, Box1DKernel
import glob

def boxcar(flux,factor):
  smoothed_flux=convolve(flux,Box1DKernel(factor))
  return smoothed_flux

def ensure_dir(d):
  if not os.path.exists(d):
      os.makedirs(d)

def filewriter(wavelength, flux, error, dq, save_path=os.getcwd()+'/', filename=''): #writes spectra to file in format (wavelength   flux   fluxerror)
  ensure_dir(save_path)
  fl=open((save_path+filename),'w') #makes a file to put the values of the spectra in
  for w, f, e, q in zip(wavelength, flux, error, dq):
	  fl.write('%f %g %g %i\n'%(w,f,e,q))

def plot_spectrum(wavelength, flux, error, dq, rootname):
	flux = boxcar(flux, 5)
	plt.figure(rootname)
	plt.plot(wavelength, flux, '0.5')
	wavelength_dq, flux_dq = wavelength[dq==0], flux[dq==0]
	plt.plot(wavelength_dq, flux_dq, label = 'spectrum')
	plt.plot(wavelength, error, label='error')  
	plt.legend()
	plt.show()
	
#stars=['cc_cet','eg_uma','lm_com','uz_sex','wd1458+171','wd1504+546','wd2317+268']
#stars=['cc_cet']
#stars = ['eggr_38']
#for i in range(len(stars)):
 #   star=stars[i]
  #  print(star)
data_loc=''
for spectrum in glob.glob(data_loc+'*x1dsum.*'):
    hdul=fits.open(spectrum)
    hdr = hdul[0].header
    hdr2 = hdul[1].header
    data = hdul[1].data
    hdul.close()

    #details of observation from headers
    #date=hdul[1].header['DATE-OBS']
   # cenwave = hdul[1].header['CENWAVE']
    #arm = hdul[1].header['DETECTOR']
    #star = hdul[1].header['TARGNAME']
    #print(arm)

#HELLO!!!
    #data
    #data=hdul[1].data
    
    w = np.array([], dtype=float)
    f = np.array([], dtype=float)
    e = np.array([], dtype=float)
    dq = np.array([], dtype=int)
    
    for j in range(len(data)):
        w = np.concatenate((w, data['wavelength'][j]))
        f = np.concatenate((f, data['flux'][j]))
        e = np.concatenate((e, data['error'][j]))
        dq = np.concatenate((dq, data['dq'][j]))

    arr1inds = w.argsort()
    w = w[arr1inds]
    f = f[arr1inds]
    e = e[arr1inds]
    dq = dq[arr1inds]

    #write to file
    filename = hdr['TARGNAME']+'_'+hdr['INSTRUME']+'_'+hdr['DETECTOR']+'_'+hdr['OPT_ELEM']+'_'+hdr2['DATE-OBS']+':'+hdr2['TIME-OBS']+'_'+hdr['ROOTNAME']+'.dat'
   # filename=(star+'_COS_'+arm+'_'+str(cenwave)+'_'+date+'.dat')
    save_path = 'spectra/'#../c25_hst/'+star+'/spectra/'  
    filewriter(w, f, e, dq, save_path = save_path, filename = filename  )

    #make a plot
    plot_spectrum(w, f, e, dq, filename)
