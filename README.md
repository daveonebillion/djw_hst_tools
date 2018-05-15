# djw_hst_tools

Repository for python scripts for working with HST spectra, based on my work for Program ID 15189. Contents so far:

* lc_extractor: extracts lightcurves from COS corrtag files

* stis_combine: stitches STIS echelle spectra together

* cos_fits_to_dat: turns COS xldsum files into ascii .dat files

## lc_extractor

lc_extractor is a script to extract a lightcurve from HST/COS corrtag files.

Tested on FUV so far, NUV to come.

Requires astropy, matplotlib, numpy.

Saves each FP_POS exposure separately, as well as one combined file.

For each exposure the counts from the A and B segments, if both present, are combined.

Airglow from Lyman alpha and ~1300A Oi is removed.

Error is photon noise only.

Optional: Plots combined lightcurve.

```python
"""
Usage: call the function lc_maker()


lc_extractor.lc_maker(star='unknown', file_path=os.getcwd()+'/', save_path=os.getcwd()+'/lightcurves/', bin_time=1., plot=True)


Arguments:
    -star = string, what you want the combined line curve to be called.
    Default is to use the 'TARGNAME' keyword in the first corrtag file it comes across.
    - file_path = string, where your corrtag files are. Default is the current directory.
    - save_path = string, where you want the output to be saves.
    Default is a new "lightcurves" directory in the current directory
    - bin_time = float, time in s to bin the lightcurve to. Default is 1s.
    - qual_check = boolean, masks out flagged pixels. Default is True.
    - plot = boolean, makes a plot of the combined lightcurve. Default is true.
    
Outputs:
    - Lightcurve of each exposure saved as [exposure rootname]_[bintime]s_lc.dat.
    - Combined lightcurve saved as [star]_[bin_time]s_lc_combined.dat.
    Lightcurves saved as time(s since MJD=0) counts(s-1) error(s-1).
"""
```
What about STIS? STIS is hard, so I may add it in the future but no promises.

## stis_combine

Stitches together STIS echelle spectra with overlapping orders.

At each overlap it interpolates the flux of the order with smaller wavelength bins onto the larger wavelength bins (although the difference is small) and makes a weighted coadd of the fluxes and errors.  

Idealy uses x1f files produced using stisblazefix (https://stisblazefix.readthedocs.io). Can use x1d files but these are badly affected by echelle blaze ripple.

For best results with E140M spectra, ask me to do a correction with the updated pht files first. This marginally improves on stisblazefix alone.

```python

"""
Usage:

call stis_echelle_coadd()

If no arguments are called, it will stitch all x1f (or x1d if no x1f files are present) in the working directory, make a new directory in the working directory and save the stiched spectra there.

Arguments:

- files: Array, list of x1f or x1d files to stitch. Default = []

- plot: boolean. Make a plot of the stitched spectrum at the end. Default = True

- nclip: int, points to clip off each echelle order to clear up the order end problems that are inherent to stis data.
         I don't know a reason to change it, but the option is there. Default =5

- file_path: string, path to where the xld/xlf files are. Remember to include a '/' at the end. Default= working directory

- save_path: string, directory to save the stitched spectra in. Remember to include a '/' at the end. Default = new directory 'stitched_spectra' in the working directory

Output:

A .dat file with space-separated columns of wavelength, flux, flux_error and data quality.

The files are named using the x1d/x1f header KEYWORDS as: 'TARGNAME_'INSTRUME_DETECTOR_OPT_ELEM_TDATEOBS:TTIMEOBS_ROOTNAME_stitched.dat'

"""
```

## cos_fits_to_dat

Extracts COS spectra from x1dsum.fits files and saves it as an ascii .dat file with columns wavelength, flux, error, dq


```python

"""
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
```