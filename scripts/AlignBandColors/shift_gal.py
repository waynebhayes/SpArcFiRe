import argparse
import glob
import multiprocessing
import random
import os
import sys
import shutil
import time
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use('agg')
from astropy.io import fits
from astropy import wcs
from psutil import virtual_memory
import PIL
from PIL import Image

import load_gals
import find_center
from galaxy import CroppingError


# Region -- Exceptions
class NoViableTemplateError(Exception): pass

class NoStarsFitError(Exception): pass

class StarsNotWithinSigmaError(Exception): pass

class FluxNotPreservedError(Exception):
    def __init__(self, input_count : float = None, output_count : float = None):
        if input_count and output_count:
            print(f"\n\nFluxNotPreservedError: Input flux count of {input_count}, output count of {output_count}.")

class ImageTooLargeError(Exception): pass

class NotEnoughViableWavebandsError(Exception):pass

class NotEnoughMemoryError(Exception):
    def __init__(self, needed_memory : float):
        print(f"\n\n{needed_memory / 1024.0 ** 3} GB of memory needed, try disabling running in parallel (-runInParallel 0) or lowering the upscale factor.\n\n")

# Used to expose exceptions when testing
class NoError(Exception): pass


# Region -- CSV Output
csv_out = None

def csv_print(*args) -> None:
    """
    Combines the args, separated by commas, and saves it to the tsv file
    """
    if csv_out is not None:
        csv_out.write(','.join([' ' if arg is None else str(arg) for arg in args]) + '\n')
    else:
        print("Unable to write to output CSV, continuing...")


# Region -- Star Processing
def filter_stars(src : "[Star]", trg : "[Star]", min_gamma : float = 1, 
                 max_gamma : float = 10, max_dist : float = np.sqrt(2)) -> "([Source Stars], [Target Stars])":
    """
    Filters out pairs of stars that are:
        -- too small (less than min_gamma)
        -- too whide (greather than max_gamma)
        -- matching star is too far (greater than max_dist)
    
    Arguments:
        src       ([Star]) : Source stars
        trg       ([Star]) : Target stars
        min_gamma (float)  : Minimum gamma (~radius) of a star to be used
        max_gamma (float)  : Maximum gamma (~radius) of a star to be used
        max_dist  (float)  : Maximum Euclidean distance between pairs of source and
                             target stars
    Returns:
        ([Star], [Star]) : (Source Stars, Target Stars) remaining after filter
    """
    m_src, m_trg = [], []
    for i in range(len(src)):
        dist = np.sqrt((src[i].x - trg[i].x)**2 + (src[i].y - trg[i].y)**2)
        if dist > max_dist: 
            continue
        m_src.append(src[i])
        m_trg.append(trg[i])

    assert len(m_src) == len(m_trg)
    return m_src, m_trg


def find_like_points(src : "[Star]", trg : "[Star]", max_dist = 5) -> "([Source Stars], [Target Stars])":
    """
    Finds stars in the target list that correspond with the source list.

    Arguments:
        src      ([Star]) : Source stars
        trg      ([Star]) : Target stars
        max_dist (float)  : Maximum Euclidean distance between the same star

    Returns:
        ([Star], [Star]) : (Source Stars, Target Stars) paired up such that Src[i] should be the same star in Trg[i]
    """
    m_src, m_trg = [], []

    for trg_star in trg:
        m_star, m_dist = None, np.inf
        for src_star in src:
            temp_dist = np.sqrt((src_star.x - trg_star.x)**2 + (src_star.y - trg_star.y)**2)
            if temp_dist < m_dist and temp_dist < max_dist:
                m_star = src_star
                m_dist = temp_dist
        
        if m_star is not None:
            m_src.append(m_star)
            m_trg.append(trg_star)

    assert len(m_src) == len(m_trg)    
    return m_src, m_trg


def fit_all_stars(galaxy : "Galaxy") -> None:
    """
    Runs curve fitting on all stars in all colorbands of the galaxy (modifies the galaxy star objects!)
    
    Arguments:
        galaxy (Galaxy) : Galaxy object with found stars present
    """
    fit_stars_dict = OrderedDict()
    for color, img, stars in galaxy.gen_img_star_pairs():
        # Calculate to subpixel accuracy the center of each star
        fit_stars = []
        for star in stars:
            try: fit_stars.append(find_center.estimate_center(img, star))
            except find_center.CurveFitError: pass
        fit_stars_dict.update({color: fit_stars})
    # Replace the galaxy stars with the fit ones
    galaxy.stars_dict = fit_stars_dict


def find_template_gal_and_stars(galaxy : "Galaxy", min_stars_template : int) -> "Color":
    """
    Determines which galaxy (if any) should be the template galaxy that all others are shifted to.
    Done by finding the galaxy with the most stars present.

    Arguments:
        galaxy             (Galaxy) : Galaxy object
        min_stars_template (int)    : Minimum needed stars present in the chosen template
    
    Returns:
        str : Color (waveband) to use
    """
    # Finds the color of the galaxy with the most stars
    template_color = sorted(galaxy.stars_dict.items(), key = lambda x: len(x[1]))[-1][0]
    # If there aren't enough stars in the template galaxy then raise error
    if len(galaxy.stars(template_color)) < min_stars_template: 
        raise NoViableTemplateError
    return template_color       



# Region -- Shifting
def average_vector(m_src : "[Star]", m_trg : "[Star]", maxsigma : float = 2) -> "[x, y]":
    """
    Returns the average vector between the source and target points, 
    excluding points outside of maxsigma * std_dev of the mean.

    Arguments:
        m_src     ([Star]) : Source stars
        m_trg     ([Star]) : Target stars
        max_sigma (float)  : Stars with a Euclidean distance of greater than max_sigma * std_dev of the
                             mean are excluded, and a new mean is calculated

    Returns:
        [x, y] : Estimated delta between the pairs of stars
    """
    # Calculate the vector between each target star and its source star
    vecs = np.array(m_src) - np.array(m_trg)
    if len(vecs) < 3:
        return np.mean(vecs, axis = 0)
    mu, dev = np.mean(vecs, axis = 0), maxsigma * np.std(vecs, axis = 0)
    # Only vectors that are within one sigma on either side of the mean
    avg = np.mean([v for v in vecs if abs(v[0] - mu[0]) < dev[0] and abs(v[1] - mu[1]) < dev[1]], axis = 0)
    return avg if not np.isnan(np.sum(avg)) else np.mean(vecs, axis = 0)


def shift_img(gal : "ndarray", vector : "[x, y]", upscale_factor : int, gal_dict : dict = dict(), color : str = "NoColor", check_count : bool = True):
    """
    Shifts the image by the given vector using Lanczos interoplation, updates the image in gal_dict

    Arguments:
        gal            (ndarray) : Image to apply the shift to
        vector         (tuple)   : Shift vector to apply to the image
        upscale_factor (int)     : Amount of times to upscale the image to for shifting
        gal_dict       (dict)    : Map of color to images to store the shifted images to
        color          (str)     : Color that is being shifted (used for storing in gal_dict)
        check_count    (bool)    : Whether or not to throw an error if the shifted image is not within 0.5% total flux count of the original

    Returns:
        ndarray : Shifted image
    """    
    print('Started processing shift {} on waveband {}...'.format(tuple(vector), color)) 
    
    if vector[0] == 0 and vector[1] == 0:
        gal_dict.update({color: gal}) 
        print('Shifted waveband {} by {} with a flux error of {} / {}'.format(color, tuple(vector), 0, np.sum(gal)))
        return gal
     
    upscale = Image.fromarray(gal)
    upscale = upscale.resize(np.array(gal.shape) * upscale_factor, resample = PIL.Image.LANCZOS)
    upscale = np.array(upscale)
    upscale = np.roll(upscale, int(round(vector[0] * upscale_factor)), axis = 1)
    upscale = np.roll(upscale, int(round(vector[1] * upscale_factor)), axis = 0)
    upscale = Image.fromarray(upscale)
    upscale = upscale.resize(gal.shape, resample = PIL.Image.LANCZOS)
    upsacle = np.array(upscale)

    # Check that the shifted output flux count is within 0.5% of the input flux count
    input_count, current_count = np.sum(gal), np.sum(upscale)
    err_perc = 0.005
    if check_count and (current_count > input_count + input_count * err_perc or current_count < input_count - input_count * err_perc):
        raise FluxNotPreservedError(input_count, current_count)
 
    print(f"Shifted waveband {color} by {tuple(vector)} with a flux error of {np.abs(input_count - current_count)} / {input_count}")
   
    gal_dict.update({color: upscale})
    return upscale


def get_galaxy_vectors(galaxy : "Galaxy", template_color : str, min_stars_all : int, ignore_sigma : bool = True) -> OrderedDict:
    """
    Calculates the shift needed for each galaxy for alignment to template_color.

    Arguments:
        galaxy         (Galaxy) : Galaxy object with images present
        template_color (str)    : Color to shift all other wavebands too match
        min_stars_all  (int)    : Minimum number of stars needed in all wavebands
        ignore_sigma   (bool)   : Whether to ignore calculated shifts of individual
                                  stars that are too far from the mean 

    Returns:
        OrderedDict : mapping of color to shift needed for alignment, shift is None if it could not be calculated
    """ 
    
    color_vectors = OrderedDict()

    for color, img, stars in galaxy.gen_img_star_pairs():
        
        # Only find matching stars if it is not the template galaxy
        if color == template_color:
            color_vectors.update({color: np.array((0, 0))})

        # Otherwise filter out stars and find the average vector between them
        else: 
            # Find stars that are near eachother in both images 
            m_src, m_trg = find_like_points(galaxy.stars_dict[template_color], stars)
            
            # If not enough points in it match the template, skip it
            if len(m_src) < min_stars_all:
                print(f"Skipping waveband {color}, could not match enough stars for realignment ({len(m_src)} stars)")
                color_vectors.update({color: None})
                continue 
            
            # Filter out stars based on thier model parameters
            m_src, m_trg = filter_stars(m_src, m_trg)

            # If less than the min stars were fit then skip this waveband
            if len(m_src) < min_stars_all:
                print(f"Skipping waveband {color}, filtered out too many stars ({len(m_src)} stars)")
                color_vectors.update({color: None})
                continue
            
            # Calculate the average vector from the reference image and this image
            try:
                color_vectors.update({color: average_vector(m_src, m_trg, ignore_sigma)})
            except StarsNotWithinSigmaError:
                print(f"Skipping waveband {color}, vector between stars disagree too much.")
                color_vectors.update({color: None})
                continue

    return color_vectors


def shift_wavebands(galaxy : "Galaxy", shift_vectors : dict, template_color : str, 
                    upscale_factor : int, run_in_parallel : bool, max_memory : int) -> OrderedDict:
    """
    Performs the shift needed for each waveband on the input images
    
    Arguments:
        galaxy          (Galaxy) : Galaxy object with input images
        shift_vecors    (dict)   : Mapping of color to needed shift for alignment
        template_color  (str)    : Color that all other wavebands are aligning to
        upscale_factor  (int)    : Amount of times each input image is upscaled for shifted
        run_in_parallel (bool)   : Whether to shift all wavebands in parallel
        max_memory      (int)    : The maximum memory (in bytes) the system is allowed to consume
    
    Returns:
        OrderedDict : Mapping of color to shifted image (stored as an ndarray)
    """
    
    shifted_imgs = multiprocessing.Manager().dict()
    procs = [multiprocessing.Process(target = shift_img, args = (galaxy.images(color), vector, upscale_factor, shifted_imgs, color)) for color, vector in shift_vectors.items() if vector is not None]

    if run_in_parallel:
        needed_memory = galaxy.width * galaxy.height * upscale_factor**2 * 8 * len(procs)
        if needed_memory > max_memory:
            raise NotEnoughMemoryError(needed_memory)
        
        for p in procs: p.start()
        for p in procs: p.join()
    
    else:
        needed_memory = galaxy.width * galaxy.height * upscale_factor**2 * 8
        print(galaxy.width, galaxy.images('g').shape)
        if needed_memory > max_memory:
            raise NotEnoughMemoryError(needed_memory)

        for p in procs:
            p.start()
            p.join()

    return shifted_imgs


# Region -- FITS Output
def save_output(outdir : str, galaxy : "Galaxy", shifted_imgs : dict, shift_vectors : dict, compressOutput: bool) -> None:
    """
    Saves the aligned galaxy images as fits files

    Arguments:
        outdir         (str)    : Output location
        galaxy         (Galaxy) : Galaxy object containing original images
        shifted_imgs   (dict)   : Map of color to shifted image (stored as ndarray)
        shift_vectors  (dict)   : Map of color to shift vector applied to each color
        compressOutput (bool)   : Whether or not to compress the output (using xz compression)
    """
    
    for color in shifted_imgs.keys():
        # Copy the header info from the original image into the new fits file
        f = galaxy.gal_dict[color]
        w = wcs.WCS(f[0].header, fix = False)
        nf = fits.PrimaryHDU()
        nf.header = f[0].header
        hdu = fits.HDUList([nf])
        nf.data = shifted_imgs[color]
        # Update shifted header  with the correct information
        try:
            hdu[0].header['CRPIX1'] += int(shift_vectors[color][0])
            hdu[0].header['CRPIX2'] += int(shift_vectors[color][1])
            hdu[0].header['CRVAL1'] += (shift_vectors[color][1] - int(shift_vectors[color][1])) * hdu[0].header['CD1_2']
            hdu[0].header['CRVAL2'] += (shift_vectors[color][0] - int(shift_vectors[color][0])) * hdu[0].header['CD2_1']
        except KeyError:
            print("Error updating header information, continuing...")

        hdu.writeto(os.path.join(outdir,  galaxy.name + '_' + color + '.fits')) 

    if compressOutput:
        for f in os.listdir(outdir):
            if not f.endswith(".xz") and f.endswith(".fits"):
                os.system(f"xz -9 -e {os.path.join(outdir, f)}")


def process_galaxy(galaxy : "Galaxy", out_dir : str, border_size : float, min_stars_template : int, min_stars_all : int, upscale_factor : int, 
                   crop_images : bool, run_in_parallel : bool, max_memory : int, compressOutput : bool, colors : tuple) -> None:
    """
    Runs AlignBandColor alignment process on the given galaxy and saves the output
    
    Arguments:
        galaxy             (Galaxy) : Input Galaxy object to align
        out_dir            (str)    : Path to save the output to
        border_size        (float)  : Percent of image size (in pixels) to make the border 
                                      (i.e. a border_size of 0.01 on an image of size 100x100 is 1 pixel wide)
        min_stars_template (int)    : Minimum number of stars needed for the template galaxy
        min_stars_all      (int)    : Minimum number of stars needed for all galaxies
        upscale_factor     (int)    : Amount to upscale the image when shifting
        crop_images        (bool)   : Whether to use Source Extractor to crop the output images down to just the galaxy
        run_in_parallel    (bool)   : Whether or not to run all image alignments in parallel
        max_memory         (int)    : Maximum memory to use (in bytes)
        compressOutput     (bool)   : Whether or not to compress the output images (using xz)
        colors             (tuple)  : The waveband colors to look for
    """
    global csv_out

    try:

        abc_path = os.path.abspath(os.path.dirname(__file__))
        tmp_path = os.path.join(abc_path, "tmp")
  
        # Set up output directory 
        p = os.path.join(out_dir, galaxy.name)
        if os.path.exists(p): 
            shutil.rmtree(p)
        os.mkdir(p)

        print("Fitting models to all stars found...")
        fit_all_stars(galaxy)
        print("Choosing reference galaxy...")
        template_color = find_template_gal_and_stars(galaxy, min_stars_template)
        template_cpy = np.copy(galaxy.images(template_color))
       
        print(f"Reference waveband chosen is {template_color} with {len(galaxy.stars_dict[template_color])} stars")
        
        shift_vectors = get_galaxy_vectors(galaxy, template_color, min_stars_all)
       
        # Crop the images to be only the galaxy
        if crop_images:
            left, right, top, bottom = galaxy.crop_images_to_galaxy()
            print(f"Cropped images down to include only the galaxy | X: ({left}, {right}) | Y: ({top}, {bottom})")
        
        # Add a border around them to preserve flux on edges
        b_size = int(galaxy.width * border_size)
        b_size = 1 if b_size < 1 else b_size
        galaxy.add_borders(b_size)
        print(f"Added border to all images of size {b_size} pixels")

        # Shift the images
        shift_imgs = shift_wavebands(galaxy, shift_vectors, template_color, upscale_factor, run_in_parallel, max_memory)    

        def len_stars(stars):
            return "NULL"  if stars is None else len(stars)

        def print_stars(stars):
            return "NULL" if stars is None else tuple((s.info() for s in stars))

        def print_vec(vec):
            return "NULL" if vec is None else str(tuple(vec))

        # Save output images and info
        save_output(p, galaxy, shift_imgs, shift_vectors, compressOutput)
        stars = galaxy.stars_dict
        
        csv_print(galaxy.name, min_stars_template, min_stars_all, upscale_factor, template_color,
                  *[print_vec(shift_vectors[c]) for c in colors],
                  *[len(stars[c])               for c in colors],
                  *[print_stars(stars[c])       for c in colors])

    except NotEnoughViableWavebandsError:
        print(f"Not enough viable wavebands ({num_viable} of {min_wavebands} needed), skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass

    except NoViableTemplateError:
        print("No viable waveband template could be found, skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass

    except FluxNotPreservedError:
        print("Error in flux preservation in shifted images, skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass
    
    except ImageTooLargeError:
        print("Input image is too large to be upscaled, try lowering the upscale factor or increasing available memory, skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass

    except CroppingError:
        print("Error cropping galaxy images, skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass

    # If any other unknown error appears, clean up then exit
    except NoError:
        print(f"Other exception encountered {e}, skipping galaxy")
        try:    shutil.rmtree(p)
        except: pass
        
    finally:
        galaxy.close()
   

# Region -- Argument Processing
if __name__ == "__main__": 

    bool_choices = ["True", "true", "T", "t", "1", "False", "false", "0", "F", "f"]

    parser = argparse.ArgumentParser()
    parser.add_argument('inDir', help = 'A directory containing the input galaxies.  Expects the directory to contain a directory for each galaxy and for each of those to contain 5 images following the name_color.fits convention.')
    
    parser.add_argument('outDir', help = 'A directory for the output images.  If it does not exist then it will be created, if it already exists then all files in it will be deleted.')
    parser.add_argument('-compressOutput', default = '1', choices = bool_choices, help = 'If true then the output fits images are compressed using xz -9 -e. Note: required xz to be available on the system.') 
    
    parser.add_argument('-starClassPerc', default = 0.7, type = float, help = 'The minimum probability confidence needed of the sextractor classification that the object is a star.  Value should be in range (0,1), default is 0.7.')
    
    parser.add_argument('-cropImages', default = '1', choices = bool_choices, help = 'If true then the input images will be cropped to only the galaxy using Source Extractor.')
    
    parser.add_argument('-borderSize', default = 0.01, type = float, help = 'Controls size of the border (as a percentage of image height in pixels) added to the image to allow room for shifting. Should not be 0 otherwise information at the edges will be lost in output images.')
    
    parser.add_argument('-minStarsTemplate', default = 3, type = int, help = 'The minimum number of stars needed in the template galaxy (the one that the other wavebands will shifted to match) for a shift to be attempted.  Default is 5')
    
    parser.add_argument('-minStarsAll', default = 1, type = int, help = 'The minimum number of stars needed in all color-bands of the galaxy for a shift to be attempted.  Any color-bands that do not satisfy this property are ignored.  Default is 2.')
    
    parser.add_argument('-upscaleFactor', default = 100, type = int, help = 'The amount that each image is upscaled using Lanczos interpolation prior to shifting.')
    
    parser.add_argument('-runInParallel', default = '0', choices = bool_choices, help = 'Will process wavebands in parallel, this requires the system to have enough memory to store all upscaled wavebands simultaneously.')

    parser.add_argument("-colorsToProcess", default = "ugriz", type = str, help = "The waveband colors to process, with each character being a color.  Defaults to 'ugriz' which correspond to wavebands u, g, r, i, and z")
    
    mem = virtual_memory().total / 1024.0**3
    parser.add_argument('-maxMemory', default = mem, type = float, help = 'The maximum amount of memory (in GB) the process can use.  At least 16GB is recommended but more will be needed for larger images and larger upscale factors.')
    
    # If not arguments provided then print out the help page
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    args = parser.parse_args() 
    
    # Check that the in directory exists and follows the format required
    if not os.path.isdir(args.inDir):
        print(args.inDir, 'is not a directory.')
        exit(1)
    
    # Check that the star_class_parc is valid
    if args.starClassPerc <= 0 or args.starClassPerc > 1:
        print('Given star class percentage', args.starClassPerc, 'does not fall within the range (0, 1].')
        exit(1)
    
    # Check that border_size is valid
    if args.borderSize < 0:
        print('Border size much be a positive value.')
        exit(1)

    # Check that upscale_factor is valid
    if args.upscaleFactor < 10:
        print('Upscale factor must be at least 10')
        exit(1)

    # Check thact the memory is valid
    if args.maxMemory <= 0:
        print('Maximum memory must be a positive value in bytes')
        exit(1)
    
    # If the output directory does not exist then create it
    try:
        if not os.path.isdir(args.outDir):
            os.mkdir(args.outDir)
    except:
        print('out_dir', args.outDir, 'could not be created.')
        exit(1)

    t = ("True", "true", "1", "T", "t")
    args.cropImages     = True if args.cropImages     in t else False
    args.runInParallel  = True if args.runInParallel  in t else False
    args.compressOutput = True if args.compressOutput in t else False
    
    colors = tuple(c for c in args.colorsToProcess)
    
    # Create output file
    csv_out = open(os.path.join(args.outDir, f"abc_out_{str(time.time()).split('.')[0]}.csv"), 'w')
    # Write header to output file
    csv_print("objID", "min_stars_template", "min_stars_all", "upscale_factor", "template_color", 
              *[f"{c}_vec"       for c in colors],
              *[f"{c}_num_stars" for c in colors],
              *[f"{c}_star_info" for c in colors])    

    
    for gal in load_gals.load_galaxies(args.inDir, args.starClassPerc, args.outDir, colors):
        # load_galaxies(...) returns the galname if it was unable to load
        if type(gal) == str: 
            print(f"Failed to load {gal}, continuing...")
            continue
        
        process_galaxy(gal, args.outDir, args.borderSize, args.minStarsTemplate, args.minStarsAll, args.upscaleFactor, args.cropImages, args.runInParallel, args.maxMemory * 1024**3, args.compressOutput, colors)
        print()

    csv_out.close()

