from astropy.io import fits
import matplotlib
matplotlib.use('agg') # needed for openlab
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
import os
import shutil
import find_center
import shift_gal
import copy
import sys
import subprocess
import load_gals
import find_center
import glob
import galaxy


class_prob = 0.9
upscale_factor = 100

def test_flux_error(outdir, gal, vectors):
    """Shifts the galaxy back and forth and checks flux and positional error"""
        
    print '\nRunning double shift tests...\n'
    
    #img = np.pad(gal[0].data[137:187, 137:187], 2, 'constant')
    img = np.pad(gal[0].data, 2, 'constant')
    outdir = os.path.abspath(outdir)
    r = np.arange(0, 1.1, 0.1)
    cycle_imgs, mean_diff, mean_diff_gal, mean_diff_bg = [], [], [], []
    org = np.copy(img)
    
    seg_img = load_gals.get_seg_img(org)
    gal_val = seg_img[int(seg_img.shape[0] / 2.0), int(seg_img.shape[1] / 2.0)]
    
    for vector in vectors:
        print 'Running delta = {}'.format(vector)
        img = shift_gal.shift_img(img, vector, upscale_factor, check_count = False)
        img = shift_gal.shift_img(img, vector * -1, upscale_factor, check_count = False)
        
        diff = np.abs(org - img)
        mean_diff.append(np.mean(diff))
        mean_diff_gal.append(np.mean(diff[seg_img == gal_val]))
        mean_diff_bg.append(np.mean(diff[seg_img == 0]))

        if vector[0] == 0.5:
            diff = img - org
            load_gals.save_fits(img, os.path.join(outdir, 'shifted_{}.fits'.format(vector)))
            load_gals.save_fits(diff, os.path.join(outdir, 'residual_{}.fits'.format(vector)))
           
            print '\n---Range in Flux of Original, Doubly-Shifted, and Residual Images---'
            
            print 'Original:                     ({}, {})'.format(np.min(org), np.max(org))
            print 'Original (Galaxy Pixels):     ({}, {})'.format(np.min(org[seg_img == gal_val]), np.max(org[seg_img == gal_val]))
            print 'Original (Background Pixels): ({}, {})\n'.format(np.min(org[seg_img == 0]), np.max(org[seg_img == 0]))
 
            
            print 'Shifted:                      ({}, {})'.format(np.min(img), np.max(img))
            print 'Shifted (Galaxy Pixels):      ({}, {})'.format(np.min(img[seg_img == gal_val]), np.max(img[seg_img == gal_val]))
            print 'Shifted (Background Pixels):  ({}, {})\n'.format(np.min(img[seg_img == 0]), np.max(img[seg_img == 0]))

            print 'Residual:                     ({}, {})'.format(np.min(diff), np.max(diff))
            print 'Residual (Galaxy Pixels):     ({}, {})'.format(np.min(diff[seg_img == gal_val]), np.max(diff[seg_img == gal_val])) 
            print 'Residual (Background Pixels): ({}, {})\n'.format(np.min(diff[seg_img == 0]), np.max(diff[seg_img == 0]))

        cycle_imgs.append(np.copy(img))
        img = np.copy(org)

    print 'Saving graphs...'

    dist = [np.sqrt(s**2 + s**2) for s in r]
    
    return mean_diff, mean_diff_gal, mean_diff_bg

    for i in range(len(cycle_imgs)):
        plt.figure()
        plt.scatter(org.flatten(), np.abs(org - cycle_imgs[i]).flatten(), s = 2)
        plt.title('Initial Flux Value vs Resulting Error ({}, {})'.format(r[i], r[i]))
        plt.xlabel('Initial Flux')
        plt.ylabel('Absolute Error in Flux')
        plt.ylim(-1, 12)
        plt.savefig(os.path.join(outdir, 'cycle_flux_pixel_error_{}.png'.format(r[i])))
        
        plt.figure()
        cutoff = 5
        plt.scatter(org[org > cutoff].flatten(), 100 * np.abs(np.abs(org - cycle_imgs[i]) / org)[org > cutoff].flatten(), s = 2)
        plt.title('Initial Flux Value vs Resulting Error ({}, {})'.format(r[i], r[i]))
        plt.xlabel('Initial Flux')
        plt.ylabel('Percent Error in Flux')
        plt.ylim(-3, 60)
        plt.savefig(os.path.join(outdir, 'cycle_flux_pixel_error_{}.png'.format(r[i])))
   
    org_stars, star_diffs = [], []
    p = os.path.join(outdir, 'temp.fits')
    load_gals.save_fits(org, p)
    
    for s in load_gals.get_sextractor_points(p):
        if s.class_prob > class_prob:
            try: org_stars.append(find_center.estimate_center(org, s))
            except: pass

    for i in range(len(r)):
        load_gals.save_fits(cycle_imgs[i], p)
        stars = []
        for s in load_gals.get_sextractor_points(p):
            if s.class_prob > class_prob:
                try: stars.append(find_center.estimate_center(cycle_imgs[i], s))
                except: pass
        src, trg = shift_gal.find_like_points(org_stars, stars)
        total_dist = sum([np.sqrt((i.x - j.x)**2 + (i.y - j.y)**2) for i, j in zip(src, trg)])
        star_diffs.append(total_dist / len(src))
    os.remove(p)

    plt.figure()
    plt.plot(dist, star_diffs)
    plt.title('Mean Difference in Star Location')
    plt.xlabel('Shift Distance in Pixels')
    plt.ylabel('Mean Location Error')
    plt.savefig(os.path.join(outdir, 'cycle_star_diff.png'), bbox_inches = 'tight')
    
    return dist, star_diffs


def test_upscale_factor(outdir, img):
     
    for vec in np.random.random((20, 2)):
        print 'Testing shift vector of {}'.format(vec)
        upscale_imgs = []
        upscales = (50, 100, 200, 300, 400)
        
        for i in upscales:
            print 'Running upscale factor {}'.format(i)
            img_cpy_str = np.copy(img)
            img_cpy_str = shift_gal.shift_img(img_cpy_str, vec, i, check_count = False)
            upscale_imgs.append(np.copy(img_cpy_str))
        
        org_stars, star_diffs = [], []
        p = os.path.join(outdir, 'temp.fits')
        load_gals.save_fits(img, p)
        for s in load_gals.get_sextractor_points(p):
            if s.class_prob > class_prob:
                try: org_stars.append(find_center.estimate_center(img, s))
                except: pass
        os.remove(p)
        for up_img in upscale_imgs:
            load_gals.save_fits(up_img, p)
            stars = []
            for s in load_gals.get_sextractor_points(p):
                if s.class_prob > class_prob:
                    try: stars.append(find_center.estimate_center(up_img, s))
                    except: pass
            
            src, trg = shift_gal.find_like_points(org_stars, stars)
            total_dist = sum([np.sqrt((s.x + vec[0] - t.x)**2 + (s.y + vec[1] - t.y)**2) for s, t in zip(src, trg)])
            star_diffs.append(total_dist / len(src))
            os.remove(p)

        plt.figure()
        plt.plot(upscales, star_diffs)
        plt.title('Mean Difference in Star Location {}'.format(vec))
        plt.xlabel('Upscale Factor')
        plt.ylabel('Mean Location Error')
        plt.savefig(os.path.join(outdir, 'upscale_single_shift_star_diff_{}.png'.format(vec)))


def test_shifts(outdir, cropped_img, img, name, random_cycles, sp_path):
    pass


def test_delta_error(outdir, gal, vectors):
    """Shifts and image by delta and then determines what delta is without knowing
       it and returns the distance between the two"""
    
    org = np.pad(np.copy(gal[0].data), 2, 'constant')
    img = np.copy(org)
    
    p = os.path.join(outdir, 'temp.fits')
    load_gals.save_fits(gal, p, isObj = True)
    org_stars = []
    for s in load_gals.get_sextractor_points(p):
        if s.class_prob > class_prob:
            try: org_stars.append(find_center.estimate_center(org, s))
            except: pass
    
    tmp_gal = load_gals.load_fits(p, returnObj = True)
    shifted_imgs = [np.copy(shift_gal.shift_img(img, vector, upscale_factor, check_count = False)) for vector in vectors]

    gal = galaxy.Galaxy({'src': gal, 'trg': tmp_gal}, {'src': org_stars, 'trg': []}, None, 'test')

    estimated_vectors = []
    for simg in shifted_imgs:
        tmp_gal[0].data = simg
        gal.gal_dict.update({'trg': tmp_gal})
        load_gals.save_fits(tmp_gal, p, isObj = True)
        stars = [s for s in load_gals.get_sextractor_points(p) if s.class_prob > class_prob]
        gal.stars_dict.update({'trg': stars})
        vecs, _ = shift_gal.get_galaxy_vectors(gal, 'src', org_stars, 2)
        estimated_vectors.append(vecs['trg'])
    return np.array([np.sqrt((v[0] + ev[0])**2 + (v[1] + ev[1])**2) for v, ev in zip(vectors, estimated_vectors)])
        

def test_positional_error(outdir, gal, vectors):
    """Shifts the image only once and checks the positional error"""
   
    print '\nRunning single shift tests...\n'
    img = np.pad(np.copy(gal[0].data), 2, 'constant')
    cycle_imgs = [] 
    org = np.copy(img)
    
    for vector in vectors:
        print 'Running delta = {}'.format(vector)
        img = shift_gal.shift_img(img, vector, upscale_factor, check_count = False)
        cycle_imgs.append(np.copy(img[2:-2, 2:-2]))
        img = np.copy(org) 
        
    dist = [np.sqrt(v[0]**2 + v[1]**2) for v in vectors]
    org_stars, star_diffs = [], []
    p = os.path.join(outdir, str(np.random.ranf()) + 'temp.fits')
    org = org[2:-2, 2:-2]
    load_gals.save_fits(gal, p, isObj = True)
     
    for s in load_gals.get_sextractor_points(p):
        if s.class_prob > class_prob:
            try: org_stars.append(find_center.estimate_center(org, s))
            except: pass

    avg_num_stars = 0
    
    try:
        for i in range(len(vectors)):
            gal[0].data = cycle_imgs[i]
            load_gals.save_fits(gal, p, isObj = True)
            stars = []
            for s in load_gals.get_sextractor_points(p):
                if s.class_prob > class_prob:
                    try: stars.append(find_center.estimate_center(cycle_imgs[i], s))
                    except: pass
            
            src, trg = shift_gal.find_like_points(org_stars, stars)
            src, trg = shift_gal.filter_stars(src, trg)
            total_dist = [np.sqrt((s.x + vectors[i][0] - t.x)**2 + (s.y + vectors[i][1]- t.y)**2) for s, t in zip(src, trg)]
            total_dist = np.sum(total_dist)
            star_diffs.append(-1 if len(src) == 0 else total_dist / len(src))
            avg_num_stars += len(src)
            print star_diffs[-1], len(src)
            if star_diffs[-1] == -1: return None, None

    except load_gals.SextractorError:
        print 'Source Extractor error, skipping....'
        return None, None

    finally:
        os.remove(p)

    return np.array(star_diffs), avg_num_stars / float(len(vectors))


def test_paper_results(single = False):
    """Testing done for the paper"""
    global upscale_factor
    
    # set up output directory
    outdir = '../paperResults'
    #if not os.path.exists(outdir): os.mkdir(outdir)
    
    # load random galaxies
    w0 = '/extra/wayne1/preserve/antholn1/SDSS_DR12'
    imgs, names  = [], []
    
    if not single:
        '''
        dirs = os.listdir(w0)
        for dir in np.random.choice(dirs, 100):
            galname = np.random.choice(os.listdir(os.path.join(w0, dir)))
            names.append(galname)
            imgs.append(load_gals.load_fits(os.path.join(w0, dir, galname, galname + '_i.fits.xz'), returnObj = True))
        '''
        #'''
        with open('names.txt') as ns:
            for n in ns.readlines()[3900:]:
                n = n.rstrip()
                p = os.path.join(w0, n[-3:], n, n + '_i.fits.xz')
                try: 
                    imgs.append(load_gals.load_fits(p, returnObj = True))
                    names.append(n)
                except: pass
                if len(imgs) == 100: break
        #'''
    else:
        pass
        #imgs, names = [load_gals.load_fits('../Other/testgal/1237648703522210023/i.fits.xz', returnObj = True)], ['test']
        #imgs, names = [load_gals.load_fits('../Other/testgal2/1237648704061636868/i.fits.xz', returnObj = True)], ['test']
   
    print 'Loaded', len(imgs), 'images for testing...'
    
    #''' 
    # test delta error within a pixel
    r = np.arange(0, 1.1, 0.1) 
    vectors = np.array([(i, j) for i in r for j in r])
    for img, name in zip(imgs, names):
        try:
            err = test_delta_error(outdir, img, vectors)
            with open(os.path.join(outdir, 'delta_err.tsv'), 'a') as f:
                f.write('\t'.join([name, '\t'.join([str(d) for d in err])]) + '\n')
            print err
        except: pass
    
    
    #'''
    '''
    # test flux error within a pixel
    r = np.arange(0, 1.01, 0.01)
    upscale_factor = 100
    vectors = np.array([i, j] for i in r for j in r])
    for img, name in zip(imgs, names):
        mean, gal, bg = test_flux_error(outdir, img, vectors)
        gal = np.array(gal)
        with open(os.path.join(outdir, 'flux_err_px.tsv'), 'a') as f:
            f.write('\t'.join([name, '\t'.join([str(d) for d in gal])]) + '\n')
  
    '''
    '''
    # test positional error within a pixel
    r = np.arange(0, 1.1, 0.1)
    upscale_factor = 10
    vectors = np.array([(i, j) for i in r for j in r])
    total_err = np.zeros(121)
    
    for img, name in zip(imgs, names):
        err, num_stars = test_positional_error(outdir, img, vectors)
        if err is not None:
            err = np.array(err)
            total_err += err
            with open(os.path.join(outdir, 'pos_err_px.tsv'), 'a') as f:
                f.write('\t'.join([name, '\t'.join([str(d) for d in err])]) + '\n')
       
    ''' 
    for f in glob.glob('core*'):
        try: os.remove(f)
        except: pass

    


