# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 13:51:57 2013

@author: darren
"""

from __future__ import division

import logging
#import matplotlib.cm as cm
#import matplotlib.pyplot as plt
import numpy as np
import os
import random
import scipy.misc
import scipy.ndimage as ndimage
import scipy.stats as stats
import shutil
import subprocess
import sys

from astropy.io import fits

logger = logging.getLogger('remove_stars_with_sextractor')

star_stdev_thres = 5.0

def remove_padding(img):
    c_rep = np.tile(np.atleast_2d(img[:, 0]).T, (1, img.shape[1]))
    r_rep = np.tile(img[0, :], (img.shape[0], 1))
#    print img.shape
#    print "c_rep shape: ", c_rep.shape
#    print "r_rep shape: ", r_rep.shape
#    pyfits.write_fits(c_rep, "c_rep.fits")
#    pyfits.write_fits(r_rep, "r_rep.fits")
    npr = np.flatnonzero(
        np.sum((img != c_rep),
            axis=1) )
    npc = np.flatnonzero(
        np.sum((img != r_rep), axis=0) )
#    print [min(npr), max(npr)]
#    print [min(npc), max(npc)]
    removed_padding = None
    orig_shape = img.shape
    if (len(npr) != img.shape[0]) or (len(npc) != img.shape[1]):
        logger.warning("Image has zero-padded rows and/or columns. Removing them.")
        r_start = np.min(npr)
        r_end = np.max(npr)+1
        c_start = np.min(npc)
        c_end = np.max(npc)+1
        img = img[r_start:r_end, c_start:c_end]
        removed_padding = ((r_start,orig_shape[0]-r_end), (c_start,orig_shape[1]-c_end))
        logger.info("the following padding was removed: {0}".format(removed_padding))
#        print img.shape
#        print np.min(npr), np.max(npr)+1
#        print np.min(npc), np.max(npc)+1
    return (img, removed_padding)
    
def restore_padding(img, pad_amts):
    # for compatibility with numpy versions before 1.7
    rpad = pad_amts[0]
    img = np.vstack((np.zeros((rpad[0], img.shape[1])), img, 
        np.zeros((rpad[1], img.shape[1]))))
    cpad = pad_amts[1]
    img = np.hstack((np.zeros((img.shape[0], cpad[0])), img, 
        np.zeros((img.shape[0], cpad[1]))))
    return img

def gen_sextractor_segmentation(in_filepath, tmp_catalog_filepath, tmp_seg_filepath):
    logger.info("running SExtractor")
    sextractor_args = ['sex', in_filepath,
	    '-c', sextractor_configfile, 
	    '-CHECKIMAGE_TYPE', 'SEGMENTATION',
	    '-CHECKIMAGE_NAME', tmp_seg_filepath,
	    '-CATALOG_NAME', tmp_catalog_filepath]
    try:
        proc = subprocess.Popen(sextractor_args, stderr=subprocess.STDOUT)
        proc.communicate()
        logger.info("SExtractor finished")
    except OSError as ose:
      logger.error("Sextractor call failed: {0}".format(ose))
      logger.info("Sextractor was called with: {0}".format(
              " ".join(sextractor_args)))
      logger.info("Check that Sextractor is available on this system.")
      return None
    except subprocess.CalledProcessError as cpe:
        logger.error("Sextractor call failed")
        logger.info("SExtractor output follows: \n{0}".format(
            cpe.output))
        logger.warning("output not produced for {0}".format(in_filename))
        return None
    logger.info("reading input and SExtractor-segmentation images")
    seg_img = fits.getdata(tmp_seg_filepath)
    return seg_img.astype(np.int)
    
def remove_nonconnected_cmpts(seg_img):
    objvals = np.sort(np.unique(seg_img))
    objvals = objvals[1:]
    for objval in objvals:
        cur_mask = (seg_img == objval)
        (cmpts, n_cmpts) = ndimage.label(cur_mask, structure=np.ones((3,3)))
        if n_cmpts > 1:
            cmpt_sizes = [np.sum(cmpts == labelnum) 
                for labelnum in range(1,n_cmpts+1)]
            keep_cmpt = np.argmax(cmpt_sizes)+1
            logger.warning(("multiple connected components found for " +
                "object {0}; keeping largest (sizes: {1})").format(
                objval, cmpt_sizes))
            seg_img[cur_mask & (cmpts != keep_cmpt)] = 0
    return seg_img
    
def clean_seg_img(seg_img):
#    plt.matshow(seg_img, cmap=cm.get_cmap('hsv')); plt.colorbar(); plt.title('before segimg cleaning')
    nobj = np.max(seg_img)
    for objnum in range(1, nobj+1):
        seg_img[ndimage.binary_fill_holes(seg_img == objnum)] = objnum
    
    weights = np.ones(9)
    weights[4] = 1.1 # break ties in favor of the center element
    def filt_fxn(arr):
        freqs = np.bincount(arr.astype(np.int), weights=weights)
        if len(freqs) > (nobj+1):
            freqs = freqs[0:-1] # don't count out-of-bounds pixels
        return np.argmax(freqs)
    
    # special value for out-of-bounds entries so they can be ignored
    padval = nobj+1 
    seg_img = ndimage.filters.generic_filter(seg_img, filt_fxn, size=3, 
        mode='constant', cval=padval)
        
#    plt.matshow(seg_img, cmap=cm.get_cmap('hsv')); plt.colorbar(); plt.title('after segimg cleaning')
#    plt.show()
    return seg_img
    
## use clean_seg_img instead
#def smooth_refobj(seg_img, refobj_num):
##    plt.matshow(seg_img, cmap=cm.ge11t_cmap('hsv')); plt.colorbar(); plt.title('before refobj smoothing')
#    strel = np.ones((3,3))
#    refobj_orig = (seg_img == refobj_num)
#    refobj_mask = ndimage.binary_opening(refobj_orig, 
#        structure=strel)
#    removed_from_refobj = refobj_orig - refobj_mask
#    
#    adj_objs = np.array([objnum for objnum in 
#        np.unique(seg_img[ndimage.binary_dilation(refobj_orig) - refobj_orig])
#        if objnum != refobj_num])
#    if len(adj_objs) > 0:
#        adj_obj_dists = np.ones(refobj_orig.shape + (len(adj_objs),))
#        for idx, objnum in enumerate(adj_objs):
#            adj_obj_dists[:, :, idx] = ndimage.distance_transform_edt(
#                seg_img != objnum)
#        closest_objs = adj_objs[np.argmin(adj_obj_dists, axis=2)]
#
#    print seg_img[removed_from_refobj]
#    seg_img[removed_from_refobj] = closest_objs[removed_from_refobj]
#    print seg_img[removed_from_refobj]
##    plt.matshow(removed_from_refobj)
##    plt.matshow(seg_img, cmap=cm.get_cmap('hsv')); plt.colorbar(); plt.title('after refobj smoothing')
##    plt.show()
#    return seg_img
    
    
def calc_adj_with_obj(seg_img, refobj_mask):
    """generate an array of adjacency proportions with the reference object.
    
    arguments:
    seg_img -- the segmentation image, with pixel values indicating object 
    membership
    refobj -- which pixels are in the reference object (boolean array with 
    same size as seg_img)
    
    returns:
    array of adjacency proportions, indexed by object number
    
    """
    
    refobj_mask = (refobj_mask > 0)
    objvals = np.unique(seg_img)
    adj_props = np.zeros(np.max(objvals) + 1)
#    plt.matshow(refobj_mask)
    for objval in objvals:
        curobj_mask = (seg_img == objval)
        surr = (ndimage.binary_dilation(ndimage.binary_fill_holes(curobj_mask)) ^ curobj_mask)
#        plt.matshow(surr);
        adj_props[objval] = np.sum(surr & refobj_mask) / np.sum(surr)
#    plt.show()
    return adj_props
    
def defragment_obj(seg_img, mainobj, min_ovl=0.0, min_brdr_stdev=-np.inf):
    logger.info("defragmenting object {0}".format(mainobj))
    removed_objects = set()
    done = False
    while not done:   
        refobj_mask = (seg_img == mainobj)
#        seg_img = smooth_refobj(seg_img, mainobj)
#        refobj_mask = ndimage.binary_opening((seg_img == mainobj), 
#            structure=np.ones((3,3)))
#        plt.matshow((seg_img == mainobj) + refobj_mask.astype(np.double))
#        plt.show()
        adj_props = calc_adj_with_obj(seg_img, refobj_mask)
        adj_props[0] = 0 # don't consider merging with background
        adj_objs = np.where(adj_props > min_ovl)[0]
        nz_adj_idxs = np.where(adj_props)
        adj_str = "(none)"
        if len(nz_adj_idxs) > 0:
            adj_str = ", ".join(
                ["{0}:{1}".format(idx, adj_props[idx]) for idx in nz_adj_idxs])
        logger.info("nonzero adjacency proportions: " + adj_str)
        done = True
        for obj in adj_objs:
            obj_mask = (seg_img == obj)
            brdr_stdev = calc_brdr_stdev(obj_mask)
#            adj_prop = adj_props[obj]
#            if adj_prop > always_incl_ovl:
#                logger.info("object {0} merged into object {1}: "
#                    "high adjacency ({2} > {3})".format(obj, mainobj,
#                    adj_prop, always_incl_ovl))
#                seg_img[obj_mask] = mainobj
#                removed_objects.add(obj)
#                done = False
            if brdr_stdev > min_brdr_stdev:
                logger.info("object {0} merged into object {1}: "
                    "sufficient adjacency ({2} > {3}) and non-circularity"
                    "({4} > {5})".format(obj, mainobj,
                    adj_props[obj], min_ovl, brdr_stdev, min_brdr_stdev))
                seg_img[obj_mask] = mainobj
                removed_objects.add(obj)
                done = False
            
        if not done:
            logger.info("object fragments found and merged; rerunning "
                "defragmentation with merged object")
        else:
            logger.info("no (new) object fragments found")
    return (seg_img, removed_objects)

def create_star_mask_rmcircular(seg_img, ctr_r, ctr_c):
    objvals = np.sort(np.unique(seg_img))
    assert objvals[0] == 0
    objvals = objvals[1:]
    logger.info("SExtractor found {0} objects".format(len(objvals)))
    
    ctrobj = seg_img[ctr_r, ctr_c]
    if ctrobj == 0:
        logger.warning('no object found at the center, looking for object nearest center')
        (rv, cv) = np.meshgrid(range(seg_img.shape[0]), range(seg_img.shape[1]))
        dists = (rv - ctr_r)**2 + (cv - ctr_c)**2
        dists[seg_img == 0] = np.inf;
        ctrobj = seg_img.flat[np.argmin[dists]]
    logger.info('segmentation image center is: ctr_r = {0}, ctr_c = {1}'.format(
        ctr_r, ctr_c))
    logger.info('center object is object {0}'.format(ctrobj))
    
    seg_img = remove_nonconnected_cmpts(seg_img)
    (seg_img, removed_objects) = defragment_obj(seg_img, ctrobj)
#    adj_props = calc_adj_with_obj(seg_img, 
#        ndimage.binary_opening(seg_img == ctrobj, structure=np.ones((3,3))))
    star_mask = np.zeros(seg_img.shape, dtype=bool)
    for objval in objvals:
        if objval in removed_objects:
            logger.info(("skipping: object {0} (it was merged into another"
                " object)").format(objval))
            continue
        cur_mask = (seg_img == objval)
#        (cmpts, n_cmpts) = ndimage.label(cur_mask, structure=np.ones((3,3)))
##        print np.unique(cmpts)
#        if n_cmpts > 1:
#            cmpt_sizes = [np.sum(cmpts == labelnum) 
#                for labelnum in range(1,n_cmpts+1)]
#            keep_cmpt = np.argmax(cmpt_sizes)+1
##            print "keep_cmpt: ", keep_cmpt
#            logger.warning(("multiple connected components found for " +
#                "object mask (sizes: {0})").format(cmpt_sizes))
#            cur_mask = (cmpts == keep_cmpt)
##            print np.count_nonzero(cur_mask)
        (all_x, all_y) = np.nonzero(cur_mask)
        mean_x = np.mean(all_x)
        mean_y = np.mean(all_y)
        mask_border = cur_mask - ndimage.binary_erosion(cur_mask)
        (brdr_x, brdr_y) = np.nonzero(mask_border)
        brdr_dists = np.sqrt((brdr_x - mean_x)**2 + (brdr_y - mean_y)**2)
        brdr_stdev = np.std(brdr_dists)
        obj_desc = ("object {0} (ctr (row, col) = ({1:.4f}, {2:.4f}), " + 
            "size = {3}, border dist stdev = {4:.4f}, " +
            "mean border dist = {5:.4f})").format(objval, mean_x, mean_y, 
            len(all_x), brdr_stdev, np.mean(brdr_dists))
        if objval == ctrobj:
            logger.info("keeping center object: " + obj_desc)
#        elif adj_props[objval] >= min_ovl:
#            logger.info("keeping object with sufficient adjacency to center: " +
#                obj_desc)
        elif brdr_stdev <= star_stdev_thres:
            star_mask |=  cur_mask
            logger.info("removing: " + obj_desc)
        else:
            logger.info("keeping: " + obj_desc)
    return star_mask

def calc_brdr_stdev(objmask):
    (all_x, all_y) = np.nonzero(objmask)
    mean_x = np.mean(all_x)
    mean_y = np.mean(all_y)
    mask_border = objmask ^ ndimage.binary_erosion(objmask)
    (brdr_x, brdr_y) = np.nonzero(mask_border)
    brdr_dists = np.sqrt((brdr_x - mean_x)**2 + (brdr_y - mean_y)**2)
    brdr_stdev = np.std(brdr_dists)
    return brdr_stdev
    
def determine_ctrobj(seg_img, ctr_r, ctr_c):
    ctrobj = seg_img[ctr_r, ctr_c]
    if ctrobj == 0:
        logger.warning('no object found at the center, looking for object nearest center')
        (rv, cv) = np.meshgrid(range(seg_img.shape[0]), range(seg_img.shape[1]))
        dists = (rv - ctr_r)**2 + (cv - ctr_c)**2
        dists[seg_img == 0] = np.inf;
        ctrobj = seg_img.flat[np.argmin(dists)]
    logger.info('segmentation image center is: ctr_r = {0}, ctr_c = {1}'.format(
        ctr_r, ctr_c))
    logger.info('center object is object {0}'.format(ctrobj))
    
    return ctrobj

def create_star_mask(seg_img, ctr_r, ctr_c):
#    dist_thres=5
    objvals = np.sort(np.unique(seg_img))
    assert objvals[0] == 0
    objvals = objvals[1:]
    logger.info("SExtractor found {0} objects".format(len(objvals)))
    
    seg_img = clean_seg_img(seg_img)
    seg_img = remove_nonconnected_cmpts(seg_img)
    ctrobj = determine_ctrobj(seg_img, ctr_r, ctr_c)
    star_mask_aggressive = (seg_img == ctrobj);
    isobj = (seg_img != 0)
    (seg_img, removed_objects) = defragment_obj(seg_img, ctrobj)
    dists_from_ctrobj = ndimage.distance_transform_edt(seg_img != ctrobj)
    star_mask = np.zeros(seg_img.shape, dtype=bool)
    for objval in objvals:
        if objval in removed_objects:
            logger.info(("skipping: object {0} (it was merged into another"
                " object)").format(objval))
            continue
        cur_mask = (seg_img == objval)
        if np.count_nonzero(cur_mask) == 0:
            logger.info(("skipping: object {0} (it has no pixels left after"
                " segmentation cleanup)").format(objval))
            continue
        min_dist = np.min(dists_from_ctrobj[cur_mask])
        brdr_stdev = calc_brdr_stdev(cur_mask)
        obj_desc = ("object {0} (dist from ctrobj = {1}, size = {2}, " +
            "border dist stdev = {3:.4f})").format(
            objval, min_dist, np.sum(cur_mask), brdr_stdev)
        if objval == ctrobj:
            logger.info("keeping center object: " + obj_desc)
        else:
            star_mask |=  cur_mask
            logger.info("removing: " + obj_desc)
#        elif min_dist <= dist_thres:
#            star_mask |=  cur_mask
#            logger.info("removing: " + obj_desc)
#        else:
#            logger.info("keeping: " + obj_desc)
    return (star_mask, star_mask_aggressive, isobj)
    
def create_star_mask_ctrcontig(seg_img, ctr_r, ctr_c):
    (ccs, n_ccs) = ndimage.label(seg_img > 0)
    ctrlbl = ccs[ctr_r, ctr_c]
    if ctrlbl == 0:
        if np.sum(ccs) == 0:
            logger.error('no components detected in image; using empty star mask')
            return np.zeros(seg_img.shape)
        logger.warning('no object at center pixel, using largest object')
        maxlbl = stats.mode(ccs[ccs != 0])[0][0]
        assert maxlbl != 0
        ctrlbl = maxlbl
        print maxlbl
    star_mask = (ccs != ctrlbl)
    return star_mask
        

# "/home/darren/Dropbox/GalaxyParser/test-data/star-removal/" test_out
# /mnt/share/in /mnt/share/out
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: remove_stars_with_sextractor in_dir out_dir"
    
    in_dirpath = sys.argv[1]
    out_dirpath = sys.argv[2]
    
    keep_seg_img = False
    write_masked_img = False
    
    logging.basicConfig(
        filename=os.path.join(out_dirpath, "remove_stars_with_sextractor.log"), 
        level=logging.DEBUG)
        
    if keep_seg_img:
        logger.info('keeping segmentation images')
    else:
        logger.info('not keeping segmentation images')
        
    if write_masked_img:
        logger.info('will write masked images (in addition to masks)')
    else:
        logger.info('will only write masks (not masked images)')
    
    if os.path.exists(out_dirpath):
        if not os.path.isdir(out_dirpath):
            errmsg = "{0} exists but is not a directory".format(out_dirpath)
            sys.stderr.write(errmsg + '\n')
            logger.fatal(errmsg)
            sys.exit(1)
        else:
            logger.warning("output directory already exists {0}".format(
                out_dirpath))
    else:
        os.mkdir(out_dirpath)
    
    sextractor_configfile = "star_rm.sex"
    if not os.path.exists(sextractor_configfile):
        errmsg = ("sextractor input files (e.g., {0}) should be in the same " + 
            "directory as this script").format(sextractor_configfile)
        sys.stderr.write(errmsg + '\n')
        logger.fatal(errmsg)
        sys.exit(1)
    
    logger.info("in directory is: {0}".format(in_dirpath))
    logger.info("out directory is: {0}".format(out_dirpath))
    
    random.seed()
    tmpdir = os.path.join(out_dirpath, "star-rm_tmp_{0:09d}".format(
        random.randrange(1, 999999999)))
    if os.path.exists(tmpdir):
        errmsg = "temporary directory {0} already exists".format(tmpdir)
        sys.stderr.write(errmsg + '\n')
        logger.error(errmsg)
        sys.exit(1)
    else:
        os.mkdir(tmpdir)
        logger.info("created temporary directory: {0}".format(tmpdir))
#    tmp_seg_filepath = os.path.join(tmpdir, "segmentation.fits")
    tmp_catalog_filepath = os.path.join(tmpdir, 'tmp_sextractor_out.cat')
    
    flist = [fname for fname in os.listdir(in_dirpath) 
        if fname.endswith(".fits")]
    logger.info("found {0} FITS files in {1}".format(len(flist), in_dirpath))
    
    fits_suffix = '.fits'
    for in_filename in flist:
        try:
            in_imgname = in_filename[:-len(fits_suffix)]
            logger.info("processing: {0}".format(in_filename))
            in_filepath = os.path.join(in_dirpath, in_filename)
            sex_in_filepath = in_filepath
            
            tmp_seg_filepath = os.path.join(
                tmpdir, in_imgname + '_segmentation.fits')
            
            in_img = fits.getdata(in_filepath)
            logger.info('input image dimensions: {0}'.format(in_img.shape))
            tmp_depad_imgpath = None
            (depad_img, removed_padding) = remove_padding(in_img)
            ctr_r = int(np.round(in_img.shape[0]/2))
            ctr_c = int(np.round(in_img.shape[1]/2))
            if removed_padding:
                logger.warning("padding found in the input image")
                tmp_depad_imgpath = os.path.join(tmpdir, 
                    in_filename[:-len(fits_suffix)] + "_depadded.fits")
                fits.writeto(tmp_depad_imgpath, depad_img)
                logger.info("created a de-padded image for SExtractor to use")
                sex_in_filepath = tmp_depad_imgpath
                ctr_r = ctr_r - removed_padding[0][0]
                ctr_c = ctr_c - removed_padding[1][0]
                ctr_r = int(round(ctr_r))
                ctr_c = int(round(ctr_c))
            logger.info('segmentation image dimensions (depadded): {0}'.format(
                depad_img.shape))
            seg_img = gen_sextractor_segmentation(sex_in_filepath, 
                tmp_catalog_filepath, tmp_seg_filepath)
            if removed_padding:
                os.remove(tmp_depad_imgpath)
                logger.info("removed temporary de-padded image")
            if not keep_seg_img:
                try:
                    os.remove(tmp_seg_filepath)
                except:
                    logger.info("not able to remove segmentation image")
                logger.info("removed segmentation image")
            if seg_img is None:
                continue
            (star_mask, star_mask_aggressive, isobj) = create_star_mask(
                seg_img, ctr_r, ctr_c)
            mask_levels = np.zeros(star_mask.shape)
            mask_levels[isobj] = 1
            mask_levels[star_mask] = 2
            mask_levels[star_mask_aggressive] = 3
            out_img = depad_img * ~star_mask
            if removed_padding:
                out_img = restore_padding(out_img, removed_padding)
                assert np.all(in_img.shape == out_img.shape)
                mask_levels = restore_padding(mask_levels, removed_padding)
            
            assert in_filename.endswith(fits_suffix)
            if write_masked_img:
                out_filepath = os.path.join(out_dirpath, in_imgname + '_star-rm.fits')
                fits.writeto(out_filepath, out_img)
                logger.info("wrote {0}".format(out_filepath))
            
            scipy.misc.imsave(os.path.join(out_dirpath, in_imgname + '_starmask.png'), mask_levels)
        except Exception as e:
            logger.warning("could not create starmask for " + in_imgname)
            logger.warning(e)
            continue
if not keep_seg_img:
    shutil.rmtree(tmpdir)
    logger.info("removed temporary directory: {0}".format(tmpdir))