from astropy.io import fits
import glob
import os
import sys
from os.path import join as pj

old_in_dir  = sys.argv[1]
old_out_dir = sys.argv[2]
new_in_dir  = sys.argv[3]

for in_fits in glob.glob(f"{old_in_dir}/*.fits"):
    gname = os.path.basename(in_fits).split(".fits")[0]
    gfits = pj(old_out_dir, gname, f"{gname}_galfit_out.fits")
    if not os.path.exists(gfits):
        continue

    with fits.open(gfits) as hdul:
        hdr  = hdul[2].header
        data = hdul[2].data
        
    fits.writeto(f"{pj(new_in_dir, os.path.basename(in_fits))}", data, hdr, overwrite=True)

