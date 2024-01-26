from astropy.io import fits
import glob
import os
import sys
import numpy as np

from os.path import join as pj

old_in_dir  = sys.argv[1]
old_out_dir = sys.argv[2]
new_in_dir  = sys.argv[3]

try:
    flip    = sys.argv[4]
except IndexError:
    print("Flip (up/down for SpArcFiRe purposes) not specified. Assuming flip.")
    flip    = "true"

for in_fits in glob.glob(f"{old_in_dir}/*.fits"):
    gname = os.path.basename(in_fits).split(".fits")[0]
    gfits = pj(old_out_dir, gname, f"{gname}_galfit_out.fits")
    if not os.path.exists(gfits):
        continue

    with fits.open(gfits) as hdul:
        hdr  = hdul[2].header
        data = hdul[2].data
        
    if flip.lower() in ("1", "true", "t"):
        data = np.flipud(data)
        
    fits.writeto(f"{pj(new_in_dir, os.path.basename(in_fits))}", data, hdr, overwrite=True)

