#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os

"""
Script to extract a lightcurve from HST/COS corrtag files
Tested on FUV so far, NUV to come.
Requires astropy, matplotlib, numpy. 
Saves each FP_POS exposure separately, as well as one combined file.
For each exposure the counts from the A and B segments, if both present, are combined.
Airglow from Lymman alpha and ~1300A Oi is removed.
Error is photon noise only. 
Optional: Plots combined lightcurve. 

Usage: call the function lc_maker()

lc_extractor.lc_maker(star='unknown', file_path=os.getcwd()+'/', save_path=os.getcwd()+'/lightcurves/', bin_time=1., plot=True)

Arguments: 
	-star = string, what you want the combined line curve to be called. 
	Default is to use the 'TARGNAME' keyword in the first corrtag file it comes across.
	- file_path = string, where your corrtag files are. Default is the curret directory.
	- save_path = sring, where you want the output to be saves. 
	Default is a new "lightcurves" directory in the current directory
	- bin_time = float, time in s to bin the lightcurve to. Default is 1.0s.
	- qual_check = boolean, masks out flagged pixels. Default is True.
	- plot = boolean, makes a plot of the combined lightcurve. Default is True.
	
Outputs: 
	- Lightcurve of each exposure saved as [exposure rootname]_[bintime]s_lc.dat.
	- Combined lightcurve saved as [star]_[bin_time]s_lc_combined.dat.
	Lightcurves saved as time(s since MJD=0) counts(s-1) error(s-1). 

"""

def region_mask(x, y, slope, intercept, height):
	mask = (y > slope*x+intercept-height/2.) & (y < slope*x+intercept+height/2.)
	return mask 
	
def ensure_dir(d):
	if not os.path.exists(d):
		os.makedirs(d)
	
def filewriter(time, counts, error, save_path, filename): 
	# writes lightcurves to dat files
	ensure_dir(save_path)
	fl=open((save_path+filename),'w')
	for t, c, e in zip(time, counts, error):
	  fl.write('%f %f %f\n'%(t, c, e))

def lc_maker(star='unknown', file_path=os.getcwd()+'/', 
             save_path=os.getcwd()+'/lightcurves/', bin_time=1., 
             qual_check=True, plot=True):

	#find the corrtag files, and end the script if there aren't any
	tag_files = glob.glob(file_path+'*corrtag*')
	if len(tag_files) == 0:
		print ('There are no corrtag files in file_path :(.')
		os._exit(1)
	
	#find all rootnames
	rootnames = np.array([], dtype=str)
	for tag in tag_files:
		rootnames= np.append(rootnames, fits.open(tag)[0].header['ROOTNAME'])
	rootnames = np.unique(rootnames)
	
	#make arrays to store combined lightcurve in
	all_time = np.array([], dtype=float)
	all_counts = np.array([], dtype=float)
	all_error = np.array([], dtype =float)
	
	for rootname in rootnames:
		
		#checks if both segments are available 
		segs = ['a', 'b']
		if (file_path+rootname+'_corrtag_a.fits') not in tag_files:
			segs = ['b']
		if (file_path+rootname+'_corrtag_b.fits') not in tag_files:
			segs = ['a']
	
		for seg in segs:
			tag_file = rootname+'_corrtag_'+seg+'.fits'
	
			seg = seg.upper() #header keywords are uppercase
			
			hdul = fits.open(file_path+tag_file)
			header = hdul[1].header
			data = hdul[1].data
			
			#get target name
			if star == 'unknown':
				star = hdul[0].header['TARGNAME']
			
			#binning to achive bin_time
			bins = int(header['EXPTIME']/bin_time)
			
			#values for extraction regions
			slope = header['SP_SLP_'+seg]
			sp_intercept = header['SP_LOC_'+seg]
			sp_height = float(header['SP_HGT_'+seg])
			
			#background regions
			bk1_intercept = header['B_BKG1_'+seg]
			bk1_height = float(header['B_HGT1_'+seg])
			bk2_intercept = header['B_BKG1_'+seg]
			bk2_height = float(header['B_HGT1_'+seg])
			
			#data 
			x = data['XCORR']
			y = data['YCORR']
			time = data['TIME']
			w = data['WAVELENGTH']
			dq = data['DQ']
			
			#mask out flagged pixels
			if qual_check == True:
				x, y, time, w = x[dq==0], y[dq==0], time[dq==0], w[dq==0]
			
			#mask out airglow from lyman alpha and oi
			wave_mask = (w < 1214.)|(w > 1217.)&(w < 1301.)|(w > 1307.)
			x, y, time = x[wave_mask], y[wave_mask], time[wave_mask]
			
			#extract lightcurve from spectrum
			sp_mask = region_mask(x, y, slope, sp_intercept, sp_height)
			sp_lc = np.histogram(time[sp_mask], bins)
			t,sp_counts = sp_lc[1][:-1], sp_lc[0]
			
			#background
			bk1_mask = region_mask(x, y, slope, bk1_intercept, bk1_height)
			bk1_lc = np.histogram(time[bk1_mask], bins)
			bk2_mask = region_mask(x, y, slope, bk2_intercept, bk2_height)
			bk2_lc = np.histogram(time[bk2_mask], bins)
			bk_counts = (bk1_lc[0]+bk2_lc[0])*(sp_height/(bk1_height+bk2_height)) #sum background counts and normalise to spectrum area
			
			#background subtraction 
			counts_bksub = sp_counts - bk_counts 
			
			#combine a and b segments, if both present
			if len(segs) > 1: 
				if seg == 'A':
					counts = counts_bksub
				else:
					counts += counts_bksub
			else:
				counts = counts_bksub
			
			#calculate photon noise
			error = counts**0.5 
			
			#convert time to absolute time
			t_adg = t + (header['EXPSTART']*86400.)
			
			#convert to counts s-1 
			counts_sec = counts/bin_time
			error_sec = error/bin_time
			
		filewriter(t_adg, counts_sec, error_sec, save_path, rootname+'_'+str(bin_time)+'s_lc.dat')
			
		all_time = np.concatenate((all_time, t_adg), axis =0)	
		all_counts = np.concatenate((all_counts, counts_sec))
		all_error = np.concatenate((all_error, error_sec))
	
	filewriter(all_time, all_counts, all_error, save_path, star+'_'+str(bin_time)+'s_lc_combined.dat') 
		
	if plot == True:
		plot_lc(star, all_time, all_counts, all_error, bin_time)

def plot_lc(star, time, counts, error, bin_time):
	plt.figure(star+'_'+str(bin_time)+'s')
	plt.subplots_adjust(top=0.99, right =0.99)
	plt.errorbar(time-time[0], counts, yerr = error, ls='none', marker='o')
	plt.xlabel('Time (s)', size=20)
	plt.ylabel('Counts (s$^{-1}$)', size=20)
	plt.show()
	

