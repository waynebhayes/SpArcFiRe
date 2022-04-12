import os
import time
import subprocess

import numpy as np
from astropy.io import fits

from galaxy import *


# Exception that is raised if an issue is encountered with Sourece Extractor
class SextractorError(Exception): pass


def get_sextractor_points(fits_path : str) -> "[Star]":
    """
    Runs the Source Extractor on the given FITS image

    Arguments:
        fits_path (str) : Path to FITS file

    Returns:
        [Star] : List of star objects found in the image
    """
    f = None 
    try:
        abc_path    = os.path.abspath(os.path.dirname(__file__))
        source_path = os.path.join(abc_path, "sex")
        tmp_path    = os.path.join(abc_path, "tmp")
        txt_path    = os.path.join(tmp_path, f"{os.path.basename(fits_path)}_star_out_{time.time()}.txt")
        
        # In case this is called from not the root directory, use absolute paths
        prev_cwd = os.getcwd()
        os.chdir(abc_path)
        proc = subprocess.Popen([source_path, fits_path, "-CATALOG_NAME", txt_path], stderr = subprocess.PIPE)
        os.chdir(prev_cwd)
        
        if proc.wait() != 0: 
            raise Exception

        stars = []
        f = open(txt_path, 'r') 
        for line in f.readlines()[4:]:
            values = ' '.join(line.rstrip().split()).split()
            stars.append(Star(float(values[0]), float(values[1]), float(values[3])))
        return stars
    
    except:
        raise SextractorError
    
    finally:
        if f is not None:
            f.close()
            try:    os.remove(txt_path)
            except: pass


def get_seg_img(img : "ndarray") -> "ndarray":
    """
    Runs Source Extractor on the given FITS image to get a segmenation image
    
    Arguments:
        img (ndarray) : 2D array of image

    Returns:
        ndarray : 2D array of segmentation image 
    """
    abc_path    = os.path.abspath(os.path.dirname(__file__))
    source_path = os.path.join(abc_path, "sex")
    tmp_path    = os.path.join(abc_path, "tmp")
        
    fits_path   = os.path.join(tmp_path, f"tmp_{time.time()}_seg.fits")
    seg_path    = os.path.join(tmp_path, f"tmp_{time.time()}_outseg.fits")
    txt_path    = os.path.join(tmp_path, f"tmp_{time.time()}_segtxt.txt")

    seg = None
    try:
        save_fits(img, fits_path)
        prev_cwd = os.getcwd()
        os.chdir(abc_path)
        proc = subprocess.Popen([source_path, fits_path, "-CHECKIMAGE_TYPE", "SEGMENTATION", "-CHECKIMAGE_NAME", seg_path, "-CATALOG_NAME", txt_path], stderr = subprocess.PIPE) 
        os.chdir(prev_cwd)

        if proc.wait() != 0: 
            raise SextractorError
    
        seg = fits.open(seg_path, ignore_missing_end = True)
        seg_img = seg[0].data
    
    except:
        raise SextractorError
    
    finally:
        if os.path.exists(fits_path) : os.remove(fits_path)
        if os.path.exists(seg_path)  : os.remove(seg_path)
        if os.path.exists(txt_path)  : os.remove(txt_path)
        if seg is not None           : seg.close()

    return seg_img


def load_fits(path : str, return_obj : bool = False) -> "ndarray or astropy.fits":
    """
    Loads the path as a fits object
    
    Arguments:
        path       (str)  : Path to load
        return_obj (bool) : Whether to return an astropy.fits object or an ndarray

    Returns:
        ndarray or astropy.fits object
    """
    abc_path = os.path.abspath(os.path.dirname(__file__))
    tmp_path = os.path.join(abc_path, "tmp")
    ext_path = os.path.join(tmp_path, f"tmp_{time.time()}_load.fits")

    # If the file is compressed, decompress it first
    if path.endswith(".xz"):    
        os.system(f"xz -d -c {path} > {ext_path}")
        path = tmp

    f = fits.open(path, ignore_missing_end = True)
    
    if return_obj: 
        return f
    else:
        img = np.copy(f[0].data)
        f.close()
        if os.path.exists(ext_path): 
            os.remove(ext_path)
        return img


def save_fits(img, path : str, is_obj : bool = False) -> None:
    """
    Saves the given image or astropy.fits object to the given path
    
    Arguments:
        img    (ndarray or astropy.fits) : Image to save
        path   (str)                     : Path to save image to
        is_obj (bool)                    : Whether img is an astropy.fits object
    """
    if is_obj:
        img.writeto(path, overwrite = True)
    else:
        fits.HDUList([fits.PrimaryHDU(data = img)]).writeto(path, overwrite = True)
 

def load_galaxy(galpath : str, galname : str, star_class_perc : float, colors : tuple) -> "Galaxy":
    """
    Loads the galaxy loacated at galpath and returns a Galaxy object

    Argumnets:
        galpath         (str)   : Path of galaxy image
        galname         (str)   : Name of the galaxy
        star_class_perc (float) : Minimum star class probability in order to be counted as a star
        colors          (tuple) : Colors to look for (i.e. 'u', 'g', etc...)

    Returns:
        Galaxy object
    """
    gal_dict = OrderedDict()
    for color in colors:
        p = os.path.join(galpath, color + '.fits')
        if os.path.exists(p):
            gal_dict.update({color: fits.open(p, ignore_missing_end = True)})
            continue
    
        p2 = os.path.join(galpath, galname + '_' + color + '.fits')    
        if os.path.exists(p2):
            gal_dict.update({color: fits.open(p2, ignore_missing_end = True)})
            continue

        p3 = os.path.join(galpath, galname + '-' + color + '.fits')
        if os.path.exists(p3):
            gal_dict.update({color: fits.open(p3, ignore_missing_end = True)})
            continue
        
        print(f"{galname} - unable to find waveband {color}, continuing with remaining...")

    # If no images were found
    if not gal_dict: return galpath
        
    # Find the stars in the galaxy images
    star_dict = OrderedDict()
    for color in gal_dict.keys():
        try:
            p = os.path.join(galpath, color + '.fits')
            if os.path.exists(p):
                star_dict.update({color: get_sextractor_points(p)})
                continue

            p2 = os.path.join(galpath, galname + '_' + color + '.fits')
            if os.path.exists(p2):    
                star_dict.update({color: get_sextractor_points(p2)})
                continue

            p3 = os.path.join(galpath, galname + '-' + color + '.fits')
            if os.path.exists(p3):
                star_dict.update({color: get_sextractor_points(p3)})
                continue
        except SextractorError:
            return galpath

    return Galaxy(gal_dict, star_dict, star_class_perc, galname)


def load_galaxies(in_dir : str, star_class_perc : float, outdir : str, colors : tuple) -> "Galaxy":
    """
    Generator that yields galaxy objects for each galaxy in in_dir.
    Assumes that in_dir exists and that the the directory structure follows what is expected.
    
    Arguments:
        in_dir          (str)   : Path to input directory
        star_class_perc (float) : Minimum star class probability in order to be counted as a star
        out_dir         (str)   : Path to output directory  
        colors          (tuple) : Colors to look for (i.e. 'u', 'g', etc...  

    Returns:
        Galaxy object if sucessfull, the galaxy name if unsucessful
    """
    
    outnames = os.listdir(outdir)
    abc_path = os.path.abspath(os.path.dirname(__file__))
    tmp_path = os.path.join(abc_path, "tmp")

    for galname in os.listdir(in_dir):
        
        print(f"--- Processing galaxy {galname} ---")

        # Skips galaxies already processed
        if galname in outnames:
            print(f"Skipping {galname}, galaxy already exists in output directory")
            continue
        
        tmpdir = None
        name = galname
        try:
            tmpdir = os.path.join(tmp_path, f"tmp_{galname}")
            if os.path.exists(tmpdir): 
                shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)

            # If the galaxy is zipped, unzip it
            files = os.listdir(os.path.join(in_dir, galname))
            if type(files) != list or len(files) < 1: 
                yield galname
            zipped = False
            if files[0].endswith(".xz"):
                zipped = True
                for f in files:
                    os.system(f"xz -d -c {os.path.join(in_dir, galname, f)} > {os.path.join(tmpdir, f[:-3])}")
            
            if zipped: 
                yield load_galaxy(tmpdir, galname, star_class_perc, colors)
            else: 
                yield load_galaxy(os.path.join(in_dir, galname), galname, star_class_perc, colors)

        except:
            yield galname
        
        finally:
            if tmpdir is not None:
                try:    shutil.rmtree(tmpdir)
                except: pass
    
