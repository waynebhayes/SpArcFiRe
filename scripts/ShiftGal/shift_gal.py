import numpy as np
import matplotlib
matplotlib.use('agg')
from astropy.io import fits
from astropy import wcs
import argparse
import glob
import os
import shutil
import sys
import load_gals
import find_center
import random
from collections import OrderedDict
import cv2
import multiprocessing
from psutil import virtual_memory
from galaxy import CroppingError

class NoViableTemplateError(Exception): pass

class NoStarsFitError(Exception): pass

class StarsNotWithinSigmaError(Exception): pass

class FluxNotPreservedError(Exception): pass

class ImageTooLargeError(Exception): pass

class NotEnoughViableWavebandsError(Exception):pass

class NotEnoughMemoryError(Exception):
    def __init__(self, needed_memory):
        print '\n\n{} GB of memory needed, try disabling running in parallel or lowering the upscale factor.\n\n'.format(needed_memory / 1024.0**3)

# used to see errors when testing
class NoError(Exception): pass

# file used to write output to
tsv_out = None


def filter_stars(src, trg, min_gamma = 1, max_gamma = 10, max_dist = np.sqrt(2)):
    """ Filters out stars that are too small, too wide, or whose matching star is too far"""
    m_src, m_trg = [], []
    for i in range(len(src)):
        dist = np.sqrt((src[i].x - trg[i].x)**2 + (src[i].y - trg[i].y)**2)
        if dist > max_dist: continue
        
        m_src.append(src[i])
        m_trg.append(trg[i])

    assert len(m_src) == len(m_trg)
    return m_src, m_trg


def find_like_points(src, trg, max_dist = 5):
    """Finds points in trg that are within max_dist distance
       from of a point in src and choses the smallest one.  
       Returns the points found in the same order."""
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


def average_vector(m_src, m_trg, maxsigma = 2):
    """Returns the average vector between the source and target points.
       It excludes points outside of maxsigma * std_dev of the mean."""
    # calculate the vector bettween each star and its source star
    vecs = np.array(m_src) - np.array(m_trg)
    if len(vecs) < 3: return np.mean(vecs, axis = 0)
    mu, dev = np.mean(vecs, axis = 0), maxsigma * np.std(vecs, axis = 0)
    # only vectors that are within one sigma on either side of the mean
    avg = np.mean([v for v in vecs if abs(v[0] - mu[0]) < dev[0] and abs(v[1] - mu[1]) < dev[1]], axis = 0)
    # if there are no points in that range, then the stars don't agree enough for the shift to be trusted
    if np.any(np.isnan(avg)): 
        raise StarsNotWithinSigmaError
    return avg


def shift_img(gal, vector, upscale_factor, gal_dict = dict(), color = 'NoColor', check_count = True):
    """Shifts the image using Lanczos interoplation, updates the image in gal_dict"""    
    print('Started processing shift {} on waveband {}...'.format(tuple(vector), color)) 
    if vector[0] == 0 and vector[1] == 0:
        gal_dict.update({color: gal}) 
        print('Shifted waveband {} by {} with a flux error of {} / {}'.format(color, tuple(vector), 0, np.sum(gal)))
        return gal
    
    padding = 1
    # if the image is large enough then its area is greater than 2^31 pixels which causes an overflow error
    # to solve this we pad the image so that it is greater than 2^32 and it believes that it is positive
    size = gal.shape[0] * upscale_factor
    if (size**2) % (2**31) >= 0 and (size**2) % (2**32) >= 2**31:
        while (size**2) % (2**31) >= 0 and (size**2) % (2**32) >= 2**31:
            size += 10
        padding = int((size - gal.shape[0] * upscale_factor) / (2 * upscale_factor)) + 1
        if padding < 1: padding = 1
        if check_count: gal = np.pad(gal, padding, 'constant')
    
    upscale = cv2.resize(gal, dsize = tuple(np.array(gal.shape) * upscale_factor), interpolation = cv2.INTER_LANCZOS4)
    upscale = np.roll(upscale, int(round(vector[0] * upscale_factor)), axis = 1)
    upscale = np.roll(upscale, int(round(vector[1] * upscale_factor)), axis = 0)
    upscale = cv2.resize(upscale, dsize = tuple(gal.shape), interpolation = cv2.INTER_LANCZOS4)

    input_count, current_count = np.sum(gal), np.sum(upscale)
    # check to make sure that the output photon count is within 0.05% of the original
    err_perc = 0.05
    if check_count and (current_count > input_count + input_count * err_perc or current_count < input_count - input_count * err_perc):
        raise FluxNotPreservedError
 
    print('Shifted waveband {} by {} with a flux error of {} / {}'.format(color, tuple(vector), np.abs(input_count - current_count), input_count))

    if padding == 1:
        gal_dict.update({color: upscale})
    else:
        upscale = upscale[padding : -1 * padding, padding : -1 * padding]
        gal_dict.update({color: upscale})

    return upscale    


def tsv_print(*args):
    """Combines the args, separated by tabs, and saves it to the tsv file"""
    tsv_out.write('\t'.join([' ' if arg is None else str(arg) for arg in args]) + '\n')


def fit_all_stars(galaxy):
    """Fits a curve to all stars in all colorbands, modifies the galaxy stars"""
    fit_stars_dict = OrderedDict()
    for color, img, stars in galaxy.gen_img_star_pairs():
        # try to calculate to subpixel accuracy the center of each star
        fit_stars = []
        for star in stars:
            try: fit_stars.append(find_center.estimate_center(img, star))
            except find_center.CurveFitError: pass
        fit_stars_dict.update({color: fit_stars})
    # replace the galaxy stars with the fit ones
    galaxy.stars_dict = fit_stars_dict


def find_template_gal_and_stars(galaxy, min_stars_template):
    """Determines which galaxy (if any) should be the template galaxy that 
       all others are shifted to.  Returns the color chosen and the fitted stars."""
    # finds the color of the galaxy with the most stars
    template_color = sorted(galaxy.stars_dict.items(), key = lambda x: len(x[1]))[-1][0]

    # if there aren't enough stars in the template galaxy then raise error
    if len(galaxy.stars(template_color)) < min_stars_template: raise NoViableTemplateError
    
    return template_color       


def get_galaxy_vectors(galaxy, template_color, min_stars_all):
    """Returns a list of 5 vectors representing the shift needed to align each to
       the template, if None then the a vector could not be found for the galaxy.""" 
    
    color_vectors = OrderedDict()

    for color, img, stars in galaxy.gen_img_star_pairs():
        
        # only find matching stars if it is not the template galaxy
        if color == template_color:
            color_vectors.update({color: np.array((0, 0))})

        # otherwise filter out stars and find the average vector between them
        else: 
            # find stars that are near eachother in both images 
            m_src, m_trg = find_like_points(galaxy.stars_dict[template_color], stars)
            
            # if not enough points in it match the template, skip it
            if len(m_src) < min_stars_all:
                print('Skipping waveband {}, could not match enough stars for realignment ({} stars)'.format(color, len(m_src)))
                color_vectors.update({color: None})
                continue 
            
            # filter out stars based on thier model parameters
            m_src, m_trg = filter_stars(m_src, m_trg)

            # if less than the min stars were fit then skip this waveband
            if len(m_src) < min_stars_all:
                print('Skipping waveband {}, filtered out too many stars ({} stars)'.format(color, len(m_src)))
                color_vectors.update({color: None})
                continue
            
            # calculate the average vector from the reference image and this image
            try:
                color_vectors.update({color: average_vector(m_src, m_trg)})
            except StarsNotWithinSigmaError:
                print('Skipping waveband {}, vector between stars disagree too much.'.format(color))
                color_vectors.update({color: None})
                continue

    return color_vectors


def shift_wavebands(galaxy, shift_vectors, template_color, upscale_factor, run_in_parallel, max_memory):
    """Returns a dict of the shifted images (only those that a vector was found for)"""
    
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
        print galaxy.width, galaxy.images('g').shape
        if needed_memory > max_memory:
            raise NotEnoughMemoryError(needed_memory)

        for p in procs:
            p.start()
            p.join()

    return shifted_imgs


def save_output(outdir, galaxy, shifted_imgs, shift_vectors, compressOutput):
    """Saves the output as fits files"""
    
    for color in shifted_imgs.keys():
       # copy the header info from the original image into the new fits file
        f = galaxy.gal_dict[color]
        w = wcs.WCS(f[0].header, fix = False)
        nf = fits.PrimaryHDU()
        nf.header = f[0].header
        hdu = fits.HDUList([nf])
        nf.data = shifted_imgs[color]
        # update shifted header and save
        hdu[0].header['CRPIX1'] += int(shift_vectors[color][0])
        hdu[0].header['CRPIX2'] += int(shift_vectors[color][1])
        hdu[0].header['CRVAL1'] += (shift_vectors[color][1] - int(shift_vectors[color][1])) * hdu[0].header['CD1_2']
        hdu[0].header['CRVAL2'] += (shift_vectors[color][0] - int(shift_vectors[color][0])) * hdu[0].header['CD2_1']
        hdu.writeto(os.path.join(outdir,  galaxy.name + '_' + color + '.fits')) 

    if compressOutput:
        for f in os.listdir(outdir):
            if not '.xz' in f and '.fits' in f:
                os.system('xz -9 -e {}'.format(os.path.join(outdir, f)))


def process_galaxy(galaxy, out_dir, border_size, min_stars_template, min_stars_all, upscale_factor, crop_images, run_in_parallel, max_memory, compressOutput):

    try:
        # set up output directory 
        p = os.path.join(out_dir, galaxy.name)
        if os.path.exists(p): shutil.rmtree(p)
        os.mkdir(p)
        global tsv_out
        
        print('--- Processing galaxy {} ---'.format(galaxy.name))

        print 'Fitting models to all stars found...'
        fit_all_stars(galaxy)
        template_color = find_template_gal_and_stars(galaxy, min_stars_template)
        template_cpy = np.copy(galaxy.images(template_color))
       
        print('Reference waveband chosen is {} with {} stars'.format(template_color, len(galaxy.stars_dict[template_color]))) 
        
        shift_vectors = get_galaxy_vectors(galaxy, template_color, min_stars_all)
       
        # crop the images to be only the galaxy
        if crop_images:
            left, right, top, bottom = galaxy.crop_images_to_galaxy()
            print 'Cropped images down to include only the galaxy | X: ({}, {}) | Y: ({}, {})'.format(left, right, top, bottom)
        
        # add a border around them to preserve flux on edges
        b_size = int(galaxy.width * border_size)
        b_size = 1 if b_size < 1 else b_size
        galaxy.add_borders(b_size)
        print 'Added border to all images of size {}'.format(b_size)

        # save averaged image for testing purposes
        avg = np.zeros(galaxy.images('i').shape)
        for img in (galaxy.images('g'), galaxy.images('i'), galaxy.images('r'), galaxy.images('z')):
            avg += (img / 4)
        load_gals.save_fits(avg, os.path.join(p, 'average_pre_shift.fits'))

        # shift the images
        shift_imgs = shift_wavebands(galaxy, shift_vectors, template_color, upscale_factor, run_in_parallel, max_memory)    

        def len_stars(stars):
            return 'NULL'  if stars is None else len(stars)

        def print_stars(stars):
            return 'NULL' if stars is None else tuple((s.info() for s in stars))

        def print_vec(vec):
            return 'NULL' if vec is None else str(tuple(vec))

        # save shifed average for testing purposes
        avg = np.zeros(galaxy.images('i').shape)
        for img in shift_imgs.values():
            avg += (img / 4)
        load_gals.save_fits(avg, os.path.join(p, 'average_shift.fits'))

        # save output images and info
        save_output(p, galaxy, shift_imgs, shift_vectors, compressOutput)
        stars = galaxy.stars_dict
        tsv_print(str(galaxy.name), min_stars_template, min_stars_all, upscale_factor, template_color,
                  print_vec(shift_vectors['g']), print_vec(shift_vectors['i']), print_vec(shift_vectors['r']), print_vec(shift_vectors['u']), print_vec(shift_vectors['z']),
                  len_stars(stars['g']), len_stars(stars['i']), len_stars(stars['r']), len_stars(stars['u']), len_stars(stars['z']),
                  print_stars(stars['g']), print_stars(stars['i']), print_stars(stars['r']), print_stars(stars['u']), print_stars(stars['z']))

    except NotEnoughViableWavebandsError:
        print 'Not enough viable wavebands ({} of {} needed), skipping galaxy'.format(num_viable, min_wavebands)
        try: shutil.rmtree(p)
        except: pass

    except NoViableTemplateError:
        print 'No viable waveband template could be found, skipping galaxy'
        try: shutil.rmtree(p)
        except: pass

    except FluxNotPreservedError:
        print 'Error in flux preservation in shifted images, skipping galaxy'
        try: shutil.rmtree(p)
        except: pass
    
    except ImageTooLargeError:
        print 'Input image is too large to be upscaled, try lowering the upscale factor or increasing available memory, skipping galaxy'
        try: shutil.rmtree(p)
        except: pass

    except CroppingError:
        print 'Error cropping galaxy images, skipping galaxy'
        try: shutil.rmtree(p)
        except: pass

    # if any other unknown error appears, clean up then exit
    except NoError:
        print 'Other exception encountered {}, skipping galaxy'.format(e)
        try: shutil.rmtree(p)
        except: pass
        
    finally:
        galaxy.close()
        try:
            for dir in glob.glob('temp*'):
                shutil.rmtree(dir)
        except: pass
   
   
if __name__ == '__main__': 
    # get all of the arguments / option
    parser = argparse.ArgumentParser()
    parser.add_argument('inDir', help = 'A directory containing the input galaxies.  Expects the directory to contain a directory for each galaxy and for each of those to contain 5 images following the name_color.fits convention.')
    parser.add_argument('outDir', help = 'A directory for the output images.  If it does not exist then it will be created, if it already exists then all files in it will be deleted.')
    parser.add_argument('-compressOutput', default = '1', choices = ['True', 'true', '1', 'False', 'false', '0'], help = 'If true then the output fits images are compressed using xz -9 -e.') 
    parser.add_argument('-starClassPerc', default = 0.7, type = float, help = 'The minimum probability confidence needed of the sextractor classification that the object is a star.  Value should be in range (0,1), default is 0.7.')
    parser.add_argument('-cropImages', default = '1', choices = ['True', 'true', '1', 'False', 'false', '0'], help = 'If true then the input images will be cropped to only the galaxy using Source Extractor.')
    parser.add_argument('-borderSize', default = 0.01, type = float, help = 'Controls size of the border (as a percentage of image height) added to the image to allow room for shifting.')
    parser.add_argument('-minStarsTemplate', default = 3, type = int, help = 'The minimum number of stars needed in the template galaxy (the one that the other wavebands will shifted to match) for a shift to be attempted.  Default is 5')
    parser.add_argument('-minStarsAll', default = 1, type = int, help = 'The minimum number of stars needed in all color-bands of the galaxy for a shift to be attempted.  Any color-bands that do not satisfy this property are ignored.  Default is 2.')
    parser.add_argument('-upscaleFactor', default = 100, type = int, help = 'The amount that each image is upscaled using Lanczos interpolation prior to shifting.')
    parser.add_argument('-runInParallel', default = '0', choices = ['True', 'true', '1', 'False', 'false', '0'], help = 'Will process wavebands in parallel, this requires the system to have enough memory to store all upscaled wavebands simultaneously.')
    mem = virtual_memory().total / 1024.0**3
    parser.add_argument('-maxMemory', default = mem, type = float, help = 'The maximum amount of memory (in GB) the process can use.  At least 16GB is recommended but more will be needed for larger images and larger upscale factors.')
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    args = parser.parse_args() 
    
    # check that the in directory exists and follows the format required
    if not os.path.isdir(args.inDir):
        print args.inDir, 'is not a directory.'
        exit(1)
    
    # check that the star_class_parc is valid
    if args.starClassPerc <= 0 or args.starClassPerc > 1:
        print 'Given star class percentage', args.starClassPerc, 'does not fall within the range (0, 1].'
        exit(1)
    
    # check that border_size is valid
    if args.borderSize < 0:
        print 'Border size much be a positive value.'
        exit(1)

    # check that upscale_factor is valid
    if args.upscaleFactor < 10:
        print 'Upscale factor must be at least 10'
        exit(1)

    # check that the memory is valid
    if args.maxMemory <= 0:
        print 'Maximum memory must be a positive value in bytes'
        exit(1)
    
    # if the output directory does not exist then create it
    try:
        if not os.path.isdir(args.outDir):
            os.mkdir(args.outDir)
    except:
        print 'out_dir', args.outDir, 'could not be created.'
        exit(1)

    t = ('True', 'true', '1')
    args.cropImages = True if args.cropImages in t else False
    args.runInParallel = True if args.runInParallel in t else False
    args.compressOutput = True if args.compressOutput in t else False
    
    tsv_out = open(os.path.join(args.outDir, 'info.tsv'), 'w')
    tsv_print('objID', 'min_stars_template', 'min_stars_all', 'upscale_factor', 'template_color', 
    'g_vec', 'i_vec', 'r_vec', 'u_vec', 'z_vec', 
    'g_num_stars', 'i_num_stars', 'r_num_stars', 'u_num_stars', 'z_num_stars',
    'g_star_info', 'i_star_info', 'r_star_info', 'u_star_info', 'z_star_info')


    for gal in load_gals.load_galaxies(args.inDir, args.starClassPerc, args.outDir):
        
        if type(gal) == str: 
            print 'Failed to load {}'.format(gal)
            continue
        
        process_galaxy(gal, args.outDir, args.borderSize, args.minStarsTemplate, args.minStarsAll, args.upscaleFactor, args.cropImages, args.runInParallel, args.maxMemory * 1024**3, args.compressOutput)
        print ''

    tsv_out.close()

