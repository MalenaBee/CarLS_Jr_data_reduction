# CarLS Jr Data Reduction Pipeline
## By Malena Bloom
This repository hosts the data reduction code used for the Carnegie Lab Spectrograph (CarLS Jr), developed during the CASSI 2025 program in support of the Via Project. 

### Necessary Packages
- numpy
- pandas
- scipy
- astropy
- matplotlib
- specreduce
- skimage
#### Reccomended Packages
- aquarel
  
### Function Descriptions
All custom functions are hosted in _SpecReduction_ -> _SpecFunctions.py_
- make1D
  - Converts spectral image (as outputted by CarLS Jr) to 1D spectrum
  - Inputs:
    - img -> FITS file image data
    - area_norm -> True (default) is normalized by area, False if normalized by max peak height
  - Returns:
    - normalized spectrum
- cropSpec
  - Crops spectral image to 1000px tall; helps reduce processing time
  - Inputs:
    - img -> FITS file image data
  - Returns:
    - cropped image
- medianImg
  - Crops all spectra in a data set, stacks the images, then finds the median
  - Inputs:
    - img_num -> number of images in data set
    - path -> path to FITS files (in form _/path/to/image_), note: drop all numbers and _.FITS_ from file name
    - rows -> pixel number in each row of the image (_FITS_IMAGE.data.shape[0]_)
    - columns -> pixel number in each column of the image (_FITS_IMAGE.data.shape[1]_)
  - Returns:
    - median image (numpy array format)
- findLine
  - finds the vertical pixel value at which the background of a spectral image ends and the spectrum begins
  - Input:
    - img -> FITS file image data
    - peak_loc -> column index of peak, found using find_peaks
    - threshold -> 2 (Default), how much brighter the next pixel must be than the previous to be recognized as a signal
  - Returns:
    - pixel value of signal start
- find2Largest
  - finds the 2 largest peaks in a 1D spectrum
  - Input:
    - peaks -> list of peak values (heights) in spectrum
    - locations -> list of peak locations
  - Returns:
    - value of largest peak
    - location of largest peak
    - value of second largest peak
    - location of second largest peak
- straightSpec
  - corrects for the systematic tilt in spectral images. Works by identifying the top of the two brightest lines in a spectral image and connecting these points in a line, then tilting the image until the slope of this line is 0.
  - Input:
    - img -> FITS file image data
    - threshold -> 2 (Default), how much brighter the next pixel must be than the previous to be recognized as a signal
  - Returns:
    - corrected image (numpy array format)
  - __Notes:__
    - This function is imperfect, so the function automatically plots the original spectrum, with the line overlaid on it and the corrected image for troubleshooting purposes
    - If the function does not work, fix by either
      1. calling the function twice, with the output of the first call used as the input of the second
      2. manually rotating the image, then calling the function to fine-tune
- findPeaks
  - finds the peaks in a 1D spectrum
  - Input:
    - spec -> 1D spectrum
    - threshold -> 0.25 (Default), height needed for peak to be recognized 
  - Returns:
    - list of locations of all peaks (lowest -> highest)
    - list of values of all peaks (order corresponds to location)
    - list of full width of all peaks (order corresponds to location)
    - list of full width half max of all peaks (order corresponds to location)
- cleanPeaks
  - eliminates double peaks from peaks list (and related lists of peak properties)
  - Input:
    - peak_props -> peak locations, peak intensities, full width, full width half max
  - Returns:
    - cleaned lists for each peak property (locations, intensities, full widths, fwhm)
- centroid
  - uses a weighted mean to centroid peaks
  - Input:
    - spec -> 1D spectrum
    - pixels -> array of pixels (length of spectrum)
    - peak_props -> peak locations, peak intensities, full width, full width half max (works best when double peaks are removed)
  - Returns:
    - list of centroid positions
    - list of residual values
  - __Notes:__
    - automatically plots both the spectrum with centroids and peaks marked and a residuals plot

    
### Typical Data Reduction Process
1. take median for stacks of bias frames, darks, and lights for both ThAr and reference UNe
2. subtract background by subtracting
   * darks - bias = _darks_
   * lights - bias = _lights_
   * _lights_ - _darks_ = cleaned_img
3. flip image
4. straighten image
5. extract 1D spectrum
6. idetify peaks, clean peak properties, and centroid spectra
7. plot the UNe spectrum (with centroids) and overlay peaks from NIST line list (__LINK REDUCED LIST__)
8. find offset between UNe data and line list]
9. rough wavelength solution using a cubic regression for offset pixels vs. wavelengths 
10. identify triplet of interest in ThAr spectrum (__DEFINE/ADD LIST__) and plot against rough wavelength solution
11. find offset between the plotted points of interest and the triplet's wavelength in NIST databases (__ADD LIST__)
12. offset the wavelength solution and plot against spectrum for final spectrum!

### Credits
Made under the guidance of Dr. Jack Piotrowski at Carnegie Science Observatories
