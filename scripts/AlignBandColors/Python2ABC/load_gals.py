import numpy as np
import os
import subprocess
import shutil
from collections import OrderedDict
from astropy.io import fits
from galaxy import Galaxy
from galaxy import Star

name = None

class SextractorError(Exception): pass

def load_fits(path, returnObj = False):
   
    if '.xz' in path:    
        tmp = os.path.join('tmp', 'temp_{}_load.fits'.format(name))
        if os.path.exists(tmp):
            os.remove(tmp)

        os.system('xz -d -c {} > {}'.format(path, tmp))
        f = fits.open(tmp, ignore_missing_end = True)
        if returnObj: return f
        img = np.copy(f[0].data)
        f.close()
        os.remove(tmp)
    else:
        f = fits.open(path, ignore_missing_end = True)
        if returnObj: return f
        img = np.copy(f[0].data)
        f.close()
    
    return img


def save_fits(img, path, isObj = False):
    if isObj:
        img.writeto(path, overwrite = True)
    else:
        fits.HDUList([fits.PrimaryHDU(data = img)]).writeto(path, overwrite = True)


def get_seg_img(img):
    """Runs Source Extractor on the given FITS image and
       returns the segmenation image"""
    
    seg = None
    save_fits(img, os.path.join('tmp', 'temp_{}_seg.fits'.format(name)))
    try:
        abc_path = os.path.abspath(os.path.dirname(__file__))
        tmp_path = os.path.join(abc_path, "tmp")
        prev_cwd = os.getcwd()
        os.chdir("./SourceExtractor")
        proc = subprocess.Popen([os.path.join(os.pardir, os.pardir, 'sex'), os.path.join(tmp_path, 'temp_{}_seg.fits'.format(name)), '-CHECKIMAGE_TYPE', 'SEGMENTATION', '-CHECKIMAGE_NAME', os.path.join(tmp_path, 'temp_{}_outseg.fits'.format(name)), '-CATALOG_NAME', os.path.join(tmp_path, 'temp_{}_seg.txt'.format(name))], stderr = subprocess.PIPE) 
        if proc.wait() != 0: 
            raise SextractorError
        os.chdir(prev_cwd)
        seg = fits.open(os.path.join(tmp_path, 'temp_{}_outseg.fits'.format(name)), ignore_missing_end = True)
        seg_img = seg[0].data
    
    except:
        raise SextractorError
    
    finally:
        if os.path.exists(os.path.join(tmp_path, 'temp_{}_seg.fits'.format(name))): os.remove(os.path.join('tmp', 'temp_{}_seg.fits'.format(name)))
        if os.path.exists(os.path.join(tmp_path, 'temp_{}_seg.txt'.format(name))): os.remove(os.path.join('tmp', 'temp_{}_seg.txt'.format(name)))
        if os.path.exists(os.path.join(tmp_path, 'temp_{}_outseg.fits'.format(name))): os.remove(os.path.join('tmp', 'temp_{}_outseg.fits'.format(name)))
        if seg is not None: seg.close()

    return seg_img


def get_sextractor_points(path):
    """Runs the sextractor on the given FITS image, returns
       an array of star objects"""
    f = None 
    try:
        abc_path = os.path.abspath(os.path.dirname(__file__))
        tmp_path = os.path.join(abc_path, "tmp")
        txt_out = os.path.join(tmp_path, '{}_star_out.txt'.format(name))
        prev_cwd = os.getcwd()
        os.chdir("./SourceExtractor")
        proc = subprocess.Popen([os.path.join(os.pardir, os.pardir, 'sex'), path, '-CATALOG_NAME', txt_out], stderr = subprocess.PIPE)
        if proc.wait() != 0: 
            raise Exception
        os.chdir(prev_cwd)
        stars = []
    
        f = open(txt_out, 'r') 
        for line in f.readlines()[4:]:
            values = ' '.join(line.rstrip().split()).split()
            stars.append(Star(float(values[0]), float(values[1]), float(values[3])))

        return stars
    
    except IndexError:
        raise SextractorError
    
    finally:
        if f is not None:
            f.close()
            try: os.remove(txt_out)
            except: pass


def load_galaxy(galpath, galname, star_class_perc):
    """Loads the galaxy loacated at galpath and returns a Galaxy object""" 
    gal_dict = OrderedDict()
    for color in ('g', 'i', 'r', 'u', 'z'):
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

    # if no images were found
    if not gal_dict: return galpath
        
    # find the stars in the galaxies
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


def load_galaxies(in_dir, star_class_perc, outdir):
    """Generator that yields galaxy objects for each galaxy in in_dir.  This assumes
       that in_dir exists and that the the directory structure follows what is expected.
       Returns a galaxy object if sucessfull, the galaxy name if unsucessful."""
    global name
    outnames = os.listdir(outdir)
    for galname in os.listdir(in_dir):
        
        # skips galaxies already processed
        if galname in outnames:
            print 'Skipping {}, galaxy already exists in output directory'.format(galname)
            continue
        
        tmpdir = None
        name = galname
        try:
            tmpdir = os.path.join('tmp', 'temp_{}'.format(galname))
            if os.path.exists(tmpdir): shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)
            # if the galaxy is zipped, unzip it
            files = os.listdir(os.path.join(in_dir, galname))
            if type(files) != list or len(files) < 1: yield galname
            zipped = False
            if '.xz' in files[0]:
                zipped = True
                for f in files:
                    os.system('xz -d -c {} > {}'.format(os.path.join(in_dir, galname, f), os.path.join(tmpdir, f[:-3])))
            
            if zipped: yield load_galaxy(tmpdir, galname, star_class_perc)
            else: yield load_galaxy(os.path.join(in_dir, galname), galname, star_class_perc)
            
        except IndexError:
            yield galname
        
        finally:
            if tmpdir is not None:
                try: shutil.rmtree(tmpdir)
                except: pass
     
