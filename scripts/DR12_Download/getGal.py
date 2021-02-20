from astropy.io import fits
from astropy import wcs
import argparse
import os
import shutil
import subprocess
import sys
import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord


class DownloadFieldsError(Exception): pass

class SextractorError(Exception): pass

class GalaxyCenterNotValidError(Exception): pass

class NoError(Exception): pass


def download_fields(RA, DEC, out_path):
    """Attempts to download the 5 waveband field images to out_path.  Returns the
       fields downloaded as fits objects in order of giruz."""
    
    #print 'Downloading {}...'.format(os.path.basename(out_path))

    proc = subprocess.Popen(['./downloadFields.sh', str(RA), str(DEC), out_path], stdout = subprocess.PIPE, universal_newlines = True)
    proc.stdout.close()
    res = proc.wait()
    if res != 0: raise DownloadFieldsError

    fields = [] 
    # return fields as astropy objects
    for f in sorted(os.listdir(out_path)):
        fields.append(fits.open(os.path.join(out_path, f), ignore_missing_end = True))
    # reverse g and i bands so that i is searched for stars first (it tends to have more)
    fields[0], fields[1] = fields[1], fields[0]
    print 'Finished downloading {} field images.'.format(len(fields))
    return fields


def crop_fits(f, refx, refy, size, path, ra, dec):
    """Crops the image to size x size, with (refx, refy)
       being the cetner of the galaxy."""
    
    assert refx >= 0; assert refy >= 0
    # check to see if size / 2 pixels are available in all directions
    # i.e. not near the edge of the image
    size = int(size / 2)
    if refx - size < 0: 
        size = refx
    if refx + size > f[0].data.shape[1] - 1:
        size = f[0].data.shape[1] - 1 - refx
    if refy - size < 0:
        size = refy
    if refy + size > f[0].data.shape[0] - 1:
        size = f[0].data.shape[0] - 1 - refy
    
    # ensure it is divisible by two (keeps ra dec correct)
    if size % 2 != 0: size -= 1

    xmin, xmax = refx - size, refx + size
    ymin, ymax = refy - size, refy + size
    
    # create cropped fits
    nf = fits.PrimaryHDU()
    nf.data = f[0].data[ymin : ymax, xmin : xmax]
    nf.header = f[0].header
    
    # update the reference pixel information to be the center of the galaxy
    nf.header['CRPIX1'] = size
    nf.header['CRPIX2'] = size
    nf.header['CRVAL1'] = float(ra)
    nf.header['CRVAL2'] = float(dec)
    return nf


def save_fits(nf, path):
    fits.HDUList([nf]).writeto(path, overwrite = True)


def get_num_stars(nf, prob):
    """Runs the sextactor on the fits image and counts the number
       of stars that have a class probability of prob."""     
    
    save_fits(nf, 'temp.fits')
    f = None 
    try:
        proc = subprocess.Popen(['./sex', 'temp.fits', '-CATALOG_NAME', 'star_out.txt'], stderr = subprocess.PIPE)
        res = proc.wait()
        if res != 0: raise Exception

        count = 0 
        f = open('star_out.txt', 'r')
        for line in f.readlines()[4:]:
            values = ' '.join(line.rstrip().split()).split()
            # only load in the points that are likely to be stars
            if float(values[3]) >= prob: count += 1 

        return count
    
    except:
        raise SextractorError
    
    finally:
        if f is not None:
            f.close()
            os.remove('star_out.txt')
        os.remove('temp.fits')


def calc_galaxy_center(header, RA, DEC):
    """Using the fits header, the center pixel of the galaxy is calculated"""
    w = WCS(header)
    coords = SkyCoord(RA, DEC, unit = 'deg')
    center_x, center_y = skycoord_to_pixel(coords, w)
    return int(center_x), int(center_y)
       
           
def save_galaxy_centered(out_path, name, RA, DEC, star_class_prob, min_num_stars):
    """ Saves the galaxy centered and with enough stars visible for realignment"""

    fields = None
    try:
        fields = download_fields(RA, DEC, out_path)
        assert len(fields) == 5
        
        '''
        ---Header Values Used---
        X Reference Pixel            - CRPIX1
        Y Reference Pixel            - CRPIX2
        X Reference Pixel Ra         - CRVAL1
        Y Reference Pixel Dec        - CRVAL2
        Ra deg change per col pixel  - CD1_1
        Ra deg change per row pixel  - CD1_2
        Dec deg change per col pixel - CD2_1
        Dec deg change per row pixel - CD2_2
        '''
        
        # i and g are reversed since i tends to have the most stars, and it will be searched first
        color_names = ('i', 'g', 'r', 'u', 'z')
        crop_size = 50
        # store the center pixels in case of re-copping
        gal_centers = []
        result_crops = []
        for i in range(5):
            
            path = os.path.join(out_path, '{}_{}.fits'.format(name, color_names[i]))
            
            center_x, center_y = calc_galaxy_center(fields[i][0].header, RA, DEC)
            gal_centers.append((center_x, center_y))
            print 'Galaxy center found to be ({}, {})'.format(center_x, center_y)
            
            # if the first waveband, figure out what size the image needs to be
            if i == 0:
                num_stars = 0
                while num_stars < min_num_stars:
                    crop_size += 50
                    cropped_nf = crop_fits(fields[i], center_x, center_y, crop_size, path, RA, DEC)

                    print cropped_nf.data.shape
                    num_stars = get_num_stars(cropped_nf, star_class_prob)
                    print '{} stars found on image size {}'.format(num_stars, int(crop_size))
                   
                    # once the image can no longer get bigger
                    if cropped_nf.data.shape[0] < int(crop_size / 2) * 2 - 10: 
                        break
             
                result_crops.append(cropped_nf)
            
            # otherwize assume the crop size has been found and just crop the image
            else:
                result_crops.append(crop_fits(fields[i], center_x, center_y, crop_size, path, RA, DEC))

        # check that all the crops are the same size, if they arent then make them the same (chose smallest one)
        for i in range(4):
            # if difference is found the find the smallest x and y crop
            if result_crops[i].data.shape[0] != result_crops[i + 1].data.shape[0]: 
                smin = np.inf
                for c in result_crops: smin = min(smin, c.data.shape[0], c.data.shape[1])
                
                print 'Recropping images to size', smin
                for f in os.listdir(out_path): 
                    os.remove(os.path.join(out_path, f))
                for i in range(5):
                    save_fits(crop_fits(fields[i], gal_centers[i][0], gal_centers[i][1], smin, path, RA, DEC), os.path.join(out_path, '{}_{}.fits'.format(name, color_names[i])))
                
                break
        else:
            for i in range(5): save_fits(result_crops[i], os.path.join(out_path, '{}_{}.fits'.format(name, color_names[i])))
               
    except NoError:
        if fields is not None: 
            for f in fields: f.close()
        
        # keep a list of all that failed
        with open(os.path.join('..', 'download_errs.txt'), 'a+') as f:
            f.write(name + '\n')
        shutil.rmtree(out_path)

    finally:
        ''' 
        # remove original (big) frame files
        if os.path.exists(out_path):
            for f in os.listdir(out_path):
                if 'frame' in f:
                    try: os.remove(os.path.join(out_path, f))
                    except: pass
        '''

        # if core dumps were generated, remove them
        for f in os.listdir('.'):
            if 'core' in f:
                try: os.remove(f)
                except: pass
                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('RA')
    parser.add_argument('DEC')
    parser.add_argument('name', help = 'The name (id) of the galaxy')
    parser.add_argument('-out_dir', default = '.', help = 'If set then the images will be saved to outdir/name')
    parser.add_argument('-overwrite', default = 'True', choices = ['True', 'true', '1', 'False', 'false', '0'], help = 'If true then if any galaxies with the same name will be overwritten')
    parser.add_argument('-min_num_stars', default = 10, type = int, help = 'The minimum number of stars needed in at least one waveband image.')
    parser.add_argument('-star_class_prob', default = 0.65, type = float, help = 'The minimum probablity that an object detected counts as a star, should be in range [0, 1].  A galaxy will require min_num_stars at this probability.')
    args = parser.parse_args()
    
    # check that star arguments are valid
    if args.min_num_stars < 0:
        print 'min_num_stars must be >= 0'
        exit(1)

    if args.star_class_prob > 1 or args.star_class_prob < 0:
        print 'star_class_prob must be in range [0, 1]'
        exit(1)

    # check that out_dir exists
    if not os.path.exists(args.out_dir):
        print out_dir, 'is not a valid directory.'
        exit(1)
    
    out_path = os.path.join(args.out_dir, args.name)
    
    # create the output directory
    if os.path.exists(out_path):
        if args.overwrite in ('True', 'true', '1'):
            shutil.rmtree(out_path)
        else:
            print out_path, 'already exists and will not be overwritten'
            exit(1)
    os.mkdir(out_path)

    save_galaxy_centered(out_path, args.name, float(args.RA), float(args.DEC), args.star_class_prob, args.min_num_stars) 
