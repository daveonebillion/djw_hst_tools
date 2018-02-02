# hst_lightcurves

lc_extractor is a script to extract a lightcurve from HST/COS corrtag files
Tested on FUV so far, NUV to come.
Requires astropy, matplotlib, numpy. 
Saves each FP_POS exposure separately, as well as one combined file.
For each exposure the counts from the A and B segments, if both present, are combined.
Airglow from Lymman alpha and ~1300A Oi is removed.
Error is photon noise only. 
Optional: Plots combined lightcurve. 

Usage: call the function lc_maker()

```python
lc_extractor.lc_maker(star='unknown', file_path=os.getcwd()+'/', save_path=os.getcwd()+'/lightcurves/', bin_time=1., plot=True)
```

Arguments: 
	-star = string, what you want the combined line curve to be called. 
	Default is to use the 'TARGNAME' keyword in the first corrtag file it comes across.
	- file_path = string, where your corrtag files are. Default is the curret directory.
	- save_path = sring, where you want the output to be saves. 
	Default is a new "lightcurves" directory in the current directory
	- bin_time = float, time in s to bin the lightcurve to. Default is 1s.
	- plot = boolean, makes a plot of the combined lightcurve. Default is true.
	
Outputs: 
	<pre>- Lightcurve of each exposure saved as [exposure rootname]_[bintime]s_lc.dat.
	- Combined lightcurve saved as [star]_[bin_time]s_lc_combined.dat.
	Lightcurves saved as time(s since MJD=0) counts(s-1) error(s-1). <\pre>

What about STIS? STIS is hard, so I may add it in the future but no promises.
