import math
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import csv
import pandas as pd

import os
from os.path import join as pj
import shutil

from joblib import Parallel, delayed

# For debugging purposes
from IPython import get_ipython
def in_notebook():
    ip = get_ipython()
    
    if ip:
        return True
    else:
        return False

_HOME_DIR = os.path.expanduser("~")
if in_notebook():
    _SPARCFIRE_DIR = pj(_HOME_DIR, "sparcfire_matt") 
    _MODULE_DIR    = pj(_SPARCFIRE_DIR, "GalfitModule")
    #_MODULE_DIR="/home/azraz/SpArcFiRe/GalfitModule"
else:
    try:
        _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
        _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
    except KeyError:
        print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
        print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
        _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
sys.path.append(_MODULE_DIR)

from Classes.Components import *
from Classes.Containers import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *
from sparc_to_galfit_feedme_gen import arc_information, galaxy_information
    
def pitch_angle(*args):
    CDEF = 0.23
    
    #args = [float(arg) for arg in args]
    r = args[0]
    theta_out = args[1]
    r_out = args[2]
    alpha = args[3]
    
    theta_out = np.radians(theta_out) #+ sky_pa)

    A = 2*CDEF/(abs(theta_out) + CDEF) - 1.00001
    B = 2 - np.arctanh(A)
    C = 2 - B
    D = theta_out * 0.5**alpha

    p1 = 0.5*(np.tanh(B*(r/r_out - 1) + 2) + 1)
    p2 = theta_out*(0.5*(r/r_out + 1))**alpha

    theta_r = p1*p2

    #print(C + B*r/r_out)
    try:
        denom = 2*r_out*(np.cosh(C + B*r/r_out))**2
    except OverflowError:
        return None, None
    dp1 = B/denom

    denom = r_out + r
    dp2 = alpha*(D/denom)*(1 + r/r_out)**alpha

    dtheta_dr = dp1*p2 + dp2*p1
    k = dtheta_dr*2

    pa = np.degrees(np.arctan(1/(abs(r*dtheta_dr))))
  
    return pa, theta_r

def sum_fourier_modes(**kwargs):
    
    very_small = 0.0000001
    #r0 = kwargs.get("r0")
    amplitudes = kwargs.get("amplitudes")
    #n = kwargs.get("n")
    phis = kwargs.get("phis")
    
    assert len(amplitudes) == len(phis), "Not enough amplitude or phi values passed to fourier mode sum function."
    
    if "theta" in kwargs:
        theta = np.radians(kwargs.get("theta"))
    else:
        xgrid = kwargs.get("xgrid")
        ygrid = kwargs.get("ygrid")
        
        xc = kwargs.get("xc") 
        yc = kwargs.get("yc") 
        
        q = kwargs.get("q")
        theta = np.arctan((ygrid - yc)/(very_small + q*xgrid - q*xc))
        
        #final_xgrid = xgrid
        #final_ygrid = ygrid
        
    
    fsum = np.array([
                     amplitude*np.cos((m + 1)*(theta + phi))
                         for m, (amplitude, phi) in enumerate(zip(amplitudes, phis))
                   ])
    
    rgrid = 1 + np.sum(fsum, axis = 0)
    return rgrid

def modify_grid(x_min,
                y_min,
                x_max, 
                y_max,
                amplitudes,
                phis,
                q,
                grid_pts = 100):
    
    xgrid = np.linspace(x_min, x_max, grid_pts)
    ygrid = np.linspace(y_min, y_max, grid_pts)
    
    #xmatrix, ymatrix = np.meshgrid(xgrid, ygrid)
    
    # Assume perfect centering
    rgrid = sum_fourier_modes(amplitudes = amplitudes,
                             phis  = phis,
                             xc    = x_max//2,
                             yc    = y_max//2,
                             q     = q,
                             xgrid = xgrid,
                             ygrid = ygrid
                            )
    
    return rgrid

def send_to_parallel(gpath, in_dir, out_dir, num_pts = 100):
    
    gname = os.path.basename(gpath)
    try:
        fits_file = OutputFits(pj(gpath, f"{gname}_galfit_out.fits"))
        #print(gname, 'ok')
    except (AttributeError, FileNotFoundError, Exception):
        print(f"Something went wrong opening galaxy {gname}! Continuing...")
        #continue
        return None

    fit_region = fits_file.feedme.header.region_to_fit

    disk = fits_file.feedme.disk
    q = disk.axis_ratio

    arms       = fits_file.feedme.arms
    fourier    = fits_file.feedme.fourier

    amplitudes = fourier.amplitudes #[v[0] for k,v in fourier.param_values.items() if k != "skip"] 
    phis       = fourier.phase_angles #[v[1] for k,v in fourier.param_values.items() if k != "skip"]

    # if disk.axis_ratio > 0.85:
    #     good_fits.pop(gname)
    #     continue

    # pitch_angles_no_fourier = []
    # thetas_no_fourier = []

    # rvals_full  = np.linspace(0.5, arms.outer_rad, num_pts) #0.5*(crop_rad[1] - crop_rad[0]), num_pts)
    # xmax = outer_rad*np.cos(arms.cumul_rot)
    # ymax = outer_rad*np.sin(arms.cumul_rot)

    xmin = 0.5
    ymin = 0.5
    xmax = fit_region[1] - fit_region[0] #outer_rad*np.cos(arms.cumul_rot)
    ymax = fit_region[3] - fit_region[2] #outer_rad*np.sin(arms.cumul_rot)
    rvals_grid  = np.linspace(xmin, 0.5*xmax, num_pts)

#        rvals_grid = np.meshgrid(rvals_full)

    # Without Fourier Modes
    # for r in rvals_full:
    #     pa, theta_r = pitch_angle(r, arms.cumul_rot, arms.outer_rad, arms.powerlaw)
    #     pitch_angles_no_fourier.append(pa)
    #     thetas_no_fourier.append(theta_r)

    # With Fourier Modes
    rgrid = np.multiply(rvals_grid, modify_grid(
                                                xmin,
                                                ymin,
                                                xmax,
                                                ymax,
                                                amplitudes,
                                                phis,
                                                q,
                                                grid_pts = num_pts
                                               )
                   )

    ones = np.ones(np.shape(rgrid))
    with_fourier = np.array(list(map(pitch_angle, rgrid, ones*arms.cumul_rot, ones*arms.outer_rad, ones*arms.powerlaw)))

    pitch_angles = with_fourier[:, 0]
    thetas = with_fourier[:, 1] - np.radians(arms.sky_position_angle) #np.radians(disk.position_angle) +

    #bad_fits[gname] = (bad_fits[gname], plotting_pa(gname, arms.cumul_rot, arms.outer_rad, arms.powerlaw, arms.sky_position_angle, arms.outer_rad, 1, 0.25, 2, 0.25, polar_coordinates, diff_pa))
    #bad_fits[gname] = plotting_pa(gname, rvals_full, pa, theta_r, arms.sky_position_angle)
    #pa_rgrid_theta[gname] = (pitch_angles, np.dstack((rvals_grid, thetas))[0,:])


    galaxy_dict = galaxy_information(gname, pj(out_dir, gname))
    crop_rad = galaxy_dict["crop_rad"]
    scale_fact_std = 2*crop_rad/256

    # Set num_arms to 1 to just grab info from the longest arc
    arc_info = arc_information(
                               gname, 
                               pj(out_dir, gname), 
                               num_arms = 2, #galaxy_dict["est_arcs"], 
                               bulge_rad = galaxy_dict["bulge_maj_axs_len"], 
                               scale_fact_std = scale_fact_std
                              )

    inner_rad = arc_info["inner_rad"]
    outer_rad = max(arc_info["outer_rad"], arms.outer_rad)

    inner_idx = np.argmin(np.abs(rvals_grid - inner_rad))
    outer_idx = np.argmin(np.abs(rvals_grid - outer_rad)) + 1

#         idx_adjust = num_pts//20
#         adjust_inner = outer_idx - idx_adjust
#         adjust_outer = outer_idx + idx_adjust

    if outer_idx >= len(pitch_angles):
        outer_idx = len(pitch_angles)

    avg_constrained_pa = np.mean(pitch_angles[inner_idx : outer_idx])
    
    #avg_constrained_pa = np.mean(pitch_angles[adjust_inner : adjust_outer])
    
    return avg_constrained_pa
      
def main(in_dir, out_dir, num_pts = 100):
    
    gpaths = glob.glob(pj(out_dir, "123*"))
    gnames = [os.path.basename(i) for i in gpaths]
    
    avg_pitch_angles = Parallel(n_jobs = -2)(
                           delayed(send_to_parallel)(
                                                 gpath,
                                                 in_dir,
                                                 out_dir,
                                                 num_pts = num_pts
                                                )
                           for gpath in gpaths
                                            )
    
    pa_dict = dict(zip(gnames, avg_pitch_angles))
    pitch_angle_df = pd.DataFrame()
        
    return pitch_angle_df.from_dict(pa_dict, orient = "index", columns = ["avg_constrained_pitch_angle"])

if __name__ == "__main__":
    
    if len(sys.argv) == 3:
        in_dir = sys.argv[1]
        out_dir = sys.argv[2]
        
    elif len(sys.argv) == 4:
        in_dir = sys.argv[1]
        # Not used but this is the usual command line input
        tmp_dir = sys.argv[2]
        out_dir = sys.argv[3]
        
    else:
        in_dir  = "sparcfire-in"
        out_dir = "sparcfire-out"
        
    if not exists(in_dir) or not exists(out_dir):
        raise(f"Cannot find input/output directories {in_dir} {out_dir}.")
    
    avg_constrained_pa = main(in_dir, out_dir)
    
    filename = pj(out_dir, "avg_constrained_pitch-angles.pkl")
    avg_constrained_pa.to_pickle(filename)
    #pd.dump(pa_rgrid_theta, open(filename, 'wb'))
    