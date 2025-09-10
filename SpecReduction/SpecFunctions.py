#GOAL: Add in check for installed packages
import numpy as np
import pandas as pd
import scipy
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from scipy import signal
import specreduce as spec
from specreduce.tracing import FlatTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
from skimage.transform import rotate
from astropy import units as u


def make1D(image, area_norm = True):
    """
    Uses specreduce to convert spectral image to 1D spectrum
    area_norm: True (default) is normalized by area, False if normalized by max peak height
    Returns normalized spectrum
    Requires the following installations:
        import specreduce as spec
        from specreduce.tracing import FlatTrace
        from specreduce.background import Background
        from specreduce.extract import BoxcarExtract
        
    image: FITS file image data
    """
    height = image.shape[0]
    spec_center = height/2 #Creates a cross-section approximately midway through the spectrum
    spec_trace = FlatTrace(image, spec_center)
    #Creates background cross-section
    if height > 1000:
        bkg_width = 20
    elif height > 500:
        bkg_width = 10
    else:
        bkg_width = 5
    bkg_sep = int((spec_center - bkg_width)/2)
    spec_bkg = Background.two_sided(image, spec_trace, separation = bkg_sep, width = bkg_width)
    #subtract!
    sub_spec = image - spec_bkg
    #extract 1D spectrum
    spec_width = bkg_width
    extraction = BoxcarExtract(sub_spec, spec_trace, width = spec_width)
    spectrum = extraction.spectrum
    #normalize spectrum
    if area_norm:
        spectrum = spectrum.flux / (1 * u.DN)
        spec_area = np.trapz(y = spectrum)
        spec_norm = spectrum / spec_area
    elif not area_norm:
        spec_norm = spectrum.flux / np.max(spectrum.flux)
    # spec_norm = np.flipud(spec_norm)
    return spec_norm
    
def cropSpec(img):
    """
    crops spectrum to image 1000 px wide
    img: image data to crop
    """
    rows = img.shape[0]
    columns = img.shape[1]
    center_y = int(rows/2)
    center_x = int(columns/2)
    position = (center_x, center_y) 
    shape = (1000, columns) 
    cutout = Cutout2D(img, position, shape)
    image = cutout.data
    return image

def medianImg(img_num, path, rows, columns):
    """
    crops spectra, stacks images, and finds the median
    img_num: number of images in stack
    path: path to files (in form /path/to/image), the program adds the necessary numbers + .FIT
    rows: amount of pixels are in each row of the image (FITS_IMAGE.data.shape[0])
    columns: amount of pixels are are in each column of the image (FITS_IMAGE.data.shape[1])
    """
    stack = np.empty((0,1000,columns))
    for x in range(0, img_num-1):
        if 0 <= x < 10:
            globals()[f'im{x}'] = fits.open(f'{path}{0}{0}{0}{x}.FIT')
        elif 10 <= x < 100:
            globals()[f'im{x}'] = fits.open(f'{path}{0}{0}{x}.FIT')
        elif 100 <= x < 1000:
            globals()[f'im{x}'] = fits.open(f'{path}{0}{x}.FIT')

    for y in range(0, img_num-1):
        image = globals()[f'im{y}'][0]
        data = image.data
        data = cropSpec(data)
        stack = np.append(stack, [data], axis = 0)

    stack_avg = np.median(stack, axis = 0)
    return stack_avg

def findLine(image, peak_loc, threshold=2):
    """
    finds the start of the spectral image (row index)
    image: FITS image data
    peak_loc: column index of peak, found using find_peaks
    threshold: how much brighter the next pixel must be than the previous to be recognized as light (default: 2)
    """
    start = np.nan
    for row in range(image.shape[0] - 1):
        pixel = image[row, peak_loc]
        pixel_below = image[row+1, peak_loc]
        if pixel >=50:
            if pixel_below >= threshold*pixel:
                start = row
                break
    return start

def find2Largest (peaks, locations):
    """
    finds the 2 largest peaks in a 1D spectrum
    peaks: list of peak values in spectrum
    locations: list of peak locations
    """
    largest = 0
    largest_loc = 0
    second = 0
    second_loc = 0
    for x in range(peaks.size):
        if float(peaks[x]) > largest:
            second = largest
            second_loc = largest_loc
            largest = float(peaks[x])
            largest_loc = locations[x]
        elif float(peaks[x]) > second:
            second = float(peaks[x])
            second_loc = locations[x]
    return largest, largest_loc, second, second_loc

def straightSpec (spec_img, threshold = 2):
    """
    Identifies the top of the brightest lines in a spectrum (image) and connects them in a line. Tilts image until the slope between the top of the lines is 0.
    If result is better, but not all the way fixed, run a second time with first correcteed image as input
    spec_img: FITS image data for spectrum
    threshold: Identifies the contrast needed between the background and the light for a line to be identified (default = 2)
    """
    #find peaks in spectrum
    spec = make1D(spec_img, False)
    loc, props = scipy.signal.find_peaks(spec, height = 0.2)
    peak_val = np.array([])
    for p in range(loc.size):
        peak_val = np.append(peak_val, spec[loc[p]])

    #find start of spectral line for brightest peaks
    peak_large, large_loc, peak_second, second_loc = find_2_largest(peak_val, loc)
    print(find_2_largest(peak_val, loc))

    large_row = find_line(spec_img, large_loc)

    second_row = find_line(spec_img, second_loc)

    #tilt spectrum
    angle = np.degrees(np.arctan((large_row - second_row)/(large_loc - second_loc)))
    corrected_spec = rotate(spec_img, angle)
    print(angle)
    
    # plot (for troubleshooting)
    plt.figure(figsize = (12,12))
    plt.imshow(spec_img, norm = 'log', vmin = 1)
    plt.title('Original Image')
    plt.colorbar(shrink = 0.4, location = 'bottom', pad = 0.04)
    plt.axhline(y = 500)
    plt.scatter([second_loc, large_loc], [second_row, large_row])
    plt.plot((second_loc, large_loc), (second_row, large_row))
    plt.show()
    
    plt.figure(figsize = (12,12))
    plt.imshow(corrected_spec, norm = 'log', vmin = 1)
    plt.title('corrected image')
    plt.colorbar(shrink = 0.4, location = 'bottom', pad = 0.04)
    plt.axhline(y = 600, color = 'b', lw = 0.5)
    plt.axhline(y = 325, color = 'b', lw = 0.5)
    plt.grid(True, color = 'b')
    ax = plt.gca()
    ax.tick_params(axis='y', right=True, labelright=True, left=True, labelleft=True)
    plt.show()
    
    return corrected_spec

def findPeaks(spec,threshold = 0.025):
    """
    Find peaks and key peak characteristics in a spectrum
    spec: 1D spectrum of interest
    threshold: Height needed for a peak to be recognized (default: 0.025)
    returns: peak locations, peak intensities, full peak width, half peak width
    """
    peak_locs, peak_props = scipy.signal.find_peaks(spec, height = threshold)
    peak_vals = peak_props['peak_heights']
    fwhm = scipy.signal.peak_widths(spec, peak_locs, rel_height = 0.5)[0]
    fw = scipy.signal.peak_widths(spec, peak_locs, rel_height = 0.95)[0]

    return peak_locs, peak_vals, fw, fwhm

def cleanPeaks(peak_props):
    """
    Eliminates double peaks from peaks list (and related lists of peak properties)
    peak_props: peak locations, peak intensities, full width, full width half max
    returns: clean peak_props
    """
    peak_locs_clean = np.array([])
    peak_vals_clean = np.array([])
    fwhm_clean = np.array([])
    fw_clean = np.array([])
    peak_locs, peak_vals, fw, fwhm = peak_props
    peak_num = len(peak_locs)
    
    for x in range(peak_num):
        low_overlap = True
        high_overlap = True
        
        min_loc = peak_locs[x] - fw[x]
        max_loc = peak_locs[x] + fw[x]

        if x > 0:
            prev_max_loc = peak_locs[x - 1] + fw[x - 1]
            if prev_max_loc < min_loc:
                low_overlap = False

        if x < peak_num - 1:
            next_min_loc = np.floor(peak_locs[x + 1] - fw[x + 1])
            if next_min_loc > max_loc:
                high_overlap = False

        if x == 0 and not high_overlap:
            peak_locs_clean = np.append(peak_locs_clean, peak_locs[x])
            peak_vals_clean = np.append(peak_vals_clean, peak_vals[x])
            fwhm_clean = np.append(fwhm_clean, fwhm[x])
            fw_clean = np.append(fw_clean, fw[x])

        elif x == peak_num - 1 and not low_overlap:
            peak_locs_clean = np.append(peak_locs_clean, peak_locs[x])
            peak_vals_clean = np.append(peak_vals_clean, peak_vals[x])
            fwhm_clean = np.append(fwhm_clean, fwhm[x])
            fw_clean = np.append(fw_clean, fw[x])

        elif not high_overlap and not low_overlap:
            peak_locs_clean = np.append(peak_locs_clean, peak_locs[x])
            peak_vals_clean = np.append(peak_vals_clean, peak_vals[x])
            fwhm_clean = np.append(fwhm_clean, fwhm[x])
            fw_clean = np.append(fw_clean, fw[x])
        else:
            continue
           
    return peak_locs_clean, peak_vals_clean, fw_clean, fwhm_clean

def centroid(spec, pixels, peak_props):
    """
    Uses the weighted mean to center pixels, then plots calculated alongside residuals plot (difference between identified peaks and centroid)
    spec: spectrum
    pixels: array of pixels (length of spectrum)
    peak_props : peak locations, peak intensities, full width, full width half max (eliminate double peaks)
    returns: centroid positions and residuals
    """
    x_center = np.array([])
    for x in range(peak_props[0].size):
        min_loc = int(np.floor(peak_props[0][x] - peak_props[3][x]))
        max_loc = int(np.ceil(peak_props[0][x] + peak_props[3][x]))
        loc_range = pixels[min_loc:max_loc]
        i_range = spec[min_loc:max_loc]
        num = 0
        for y in range(loc_range.size):
            num += loc_range[y]*i_range[y]
           
        x_center = np.append(x_center, num/np.sum(i_range))

    fig, ax = plt.subplots(1, 2, figsize = (12, 4))
    fig.suptitle('Centroiding Peaks')
    
    ax[0].plot(spec, label = 'spectrum') #plots original data
    for x in range(x_center.size - 1):
        ax[0].axvline(x_center[x], color = '#bf616a', lw = 1) #centroid
    ax[0].axvline(x_center[x_center.size - 1], color = '#bf616a', lw = 1, label = 'centroid')
    ax[0].scatter(peak_props[0], peak_props[1], color = '#b48ead', lw = 0.5, label = 'peak')#original peaks
    ax[0].set_title('Spectrum with Peaks and Centroids')
    ax[0].legend()
    ax[0].set_xlabel('x Pixels')
    ax[0].set_ylabel('Normalized Flux')
    ax[0].set_ylim(0,1)
    
    residuals = np.array([])
    for p in range(peak_props[0].size):
        residual = peak_props[0][p] - x_center[p]
        residuals = np.append(residuals, residual)
        
    ax[1].scatter(np.arange(len(residuals)),residuals)
    ax[1].axhline(0,ls='-.', c = '#bf616a')
    ax[1].set_title('Centroiding Residuals')
    ax[1].set_xlabel('Peak Number')
    ax[1].set_ylabel('Residual (Peak Value - Centroid Value)')

    return x_center, residuals
