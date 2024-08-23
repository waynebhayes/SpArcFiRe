import math
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
#import glob
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

from Classes.Parameters import *
from Classes.Components import *
from Classes.Containers import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *
from sparc_to_galfit_feedme_gen import arc_information, galaxy_information
    
def pitch_angle(*args):
    CDEF = 0.23
    
    #args = [float(arg) for arg in args]
    r         = args[0]
    theta_out = args[1]
    r_out     = args[2]
    alpha     = args[3]
    r_in      = np.zeros(np.shape(r))
    if len(args) == 5:
        r_in  = args[4]
    
    theta_out = np.radians(theta_out) #+ sky_pa)

    A = 2*CDEF/(np.abs(theta_out) + CDEF) - 1.00001
    B = (2 - np.arctanh(A))*(r_out/(r_out - r_in))
    C = 2 - B
    D = theta_out * 0.5**alpha
    inv_r_out = 1/r_out

    p1 = 0.5*(np.tanh(B*(inv_r_out*r - 1) + 2) + 1)
    p2 = theta_out*(0.5*(inv_r_out*r + 1))**alpha

    theta_r = p1*p2

    try:
        denom = 2*r_out*(np.cosh(C + B*inv_r_out*r))**2
    except OverflowError:
        return None, None
    
    dp1 = B/denom

    denom = r_out + r
    dp2 = alpha*(D/denom)*(1 + inv_r_out*r)**alpha

    dtheta_dr = dp1*p2 + dp2*p1
    k = dtheta_dr*2

    pa = np.degrees(np.arctan2(1, abs(r*dtheta_dr)))
  
    return pa, theta_r

def sum_fourier_modes(fourier_dict, **kwargs):
    
    very_small = 0.0000001
    #r0 = kwargs.get("r0")
    #amplitudes = kwargs.get("amplitudes")
    #n = kwargs.get("n")
    #phis = kwargs.get("phis")
    
    #assert len(amplitudes) == len(phis), "Not enough amplitude or phi values passed to fourier mode sum function."
    
    if "theta" in kwargs:
        theta = np.radians(kwargs.get("theta"))
    else:
        default_0 = 0.0
        
        # Defaults give us a fourier mode function we can unit test if need be
        xgrid = kwargs.get("xgrid", np.linspace(default_0, 100, 100))
        ygrid = kwargs.get("ygrid", np.linspace(default_0, 100, 100))
        
        x0 = kwargs.get("x0", default_0) 
        y0 = kwargs.get("y0", default_0) 
        
        q = kwargs.get("q", 1)
        theta = np.arctan2(ygrid - y0, very_small + q*xgrid - q*x0)
                
    rotate = True
    
    fsum = np.zeros(np.shape(theta))
    for Fstr, Fmode in fourier_dict.items():
        if Fstr == "skip":
            continue
        
        amplitude = Fmode.amplitude
        phi       = Fmode.phase_angle
        
        m = float(Fstr[-1])
        phi = np.radians(phi)
        
        if m % 2 == 1 and rotate:
            theta += 0.5*np.pi
            rotate = False
            
        #print(f"Mode: {m}, amp: {amplitude}, angle: {phi}")
        
        fsum += amplitude*np.cos(m*(theta + phi))
    
    rgrid = 1 + fsum #np.sum(np.array(fsum), axis = 0)
    return rgrid #, theta

#     fsum = np.array([
#                      amplitude*np.cos(float(Fstr[-1])*(theta + np.radians(phi)))
#                          for Fstr, (amplitude, phi) in fourier_dict.items()
#                    ])
    
#     rgrid = 1 + np.sum(fsum, axis = 0)
#    return rgrid

# def modify_grid(x_min,
#                 y_min,
#                 x_max, 
#                 y_max,
#                 amplitudes,
#                 phis,
#                 q,
#                 grid_pts = 100):
    
#     xgrid = np.linspace(x_min, x_max, grid_pts)
#     ygrid = np.linspace(y_min, y_max, grid_pts)
    
#     #xmatrix, ymatrix = np.meshgrid(xgrid, ygrid)
    
#     # Assume perfect centering
#     rgrid = sum_fourier_modes(amplitudes = amplitudes,
#                              phis  = phis,
#                              xc    = x_max//2,
#                              yc    = y_max//2,
#                              q     = q,
#                              xgrid = xgrid,
#                              ygrid = ygrid
#                             )
    
#     return rgrid

def plot_scatter(
                 pitch_angles, 
                 rgrid, 
                 galaxy_info, 
                 arc_info, 
                 inner_idx = 0, 
                 outer_idx = -1, 
                 scatter_dir = "./"
                 ):
    
    ones = np.ones(len(rgrid))
    gname = str(galaxy_info.name[0])
    
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.scatter(rgrid, pitch_angles)#, label='_nolegend_')
    plt.axvline(x = rgrid[inner_idx], color = 'mediumseagreen', alpha = 0.5)
    plt.axvline(x = rgrid[outer_idx], color = 'mediumseagreen', alpha = 0.5, label='_nolegend_')

    sparc_pa  = float(abs(galaxy_info[' pa_alenWtd_avg_domChiralityOnly'].iloc[0]))
    sparc_unc = np.degrees(np.arctan2(float(arc_info.loc[0, 'num_pixels']),float(arc_info.loc[0, 'arc_length'])**2))

    inner_rad = rgrid[inner_idx]
    outer_rad = rgrid[outer_idx]
    
    limited_rgrid = rgrid[(rgrid >= inner_rad) & (rgrid <= outer_rad)]
    ones = np.ones(len(limited_rgrid))
    
    plt.plot(limited_rgrid, ones*sparc_pa, linewidth=1, color = 'orange', alpha = 0.7)
    
    # Error bars
    
    plt.plot(limited_rgrid, ones*(sparc_pa + sparc_unc), linestyle = "--", color = 'orange', alpha = 0.3, label='_nolegend_')
    
    plt.plot(limited_rgrid, ones*(sparc_pa - sparc_unc), linestyle = "--", color = 'orange', alpha = 0.3, label='_nolegend_')
    
    plt.fill_between(sorted(limited_rgrid), sparc_pa - sparc_unc, sparc_pa + sparc_unc, facecolor='orange', alpha=.2, label='_nolegend_')
    

    plt.legend(["Pitch Angle", "Longest Arc Inner/Outer Rad", "SpArcFiRe PA +/- Uncertainty"], loc = "upper right")

    plt.title(f'Galaxy {gname}')
    plt.xlabel(f'Radius')
    plt.ylabel('Pitch Angle (degrees)')
    
    #plt.xlim([0, ])
    #plt.ylim([0, ])
    #plt.show()

    filename = pj(scatter_dir, f"{gname}_scatter.png")
    plt.savefig(filename)
    plt.close()
    return filename

def plot_validation(pitch_angles, 
                    rgrid, 
                    thetas, 
                    galaxy_info, 
                    arc_info, 
                    fits_file_obj,
                    model_obj,
                    inner_idx = 0, 
                    outer_idx = -1, 
                    radial_steps = 5, 
                    validation_dir = "./",
                    **kwargs
                   ):
    
    filename = ""
    
    ones = np.ones(len(rgrid))
    gname = str(galaxy_info.name[0])
    num_pts = len(rgrid)
    
    fit_region = fits_file_obj.feedme.header.region_to_fit
    xmax = fit_region.x2 - fit_region.x1
    ymax = fit_region.y2 - fit_region.y1
    
    rval_max = np.max(rgrid)
    
    x0 = kwargs.get("x0", fits_file_obj.feedme.disk.position.x - fit_region.x1) #0.5*xmax
    y0 = kwargs.get("y0", fits_file_obj.feedme.disk.position.y - fit_region.y1) #0.5*ymax
    
    x_vals = rgrid*np.cos(thetas) + x0
    y_vals = rgrid*np.sin(thetas) + y0
    
    x_vals_flip = -rgrid*np.cos(thetas) + x0
    y_vals_flip = -rgrid*np.sin(thetas) + y0
    
    # TODO: No inclination for now
    #incl = np.radians(fits_file_obj.feedme.arms.inclination)
    #y_vals = y_vals*np.cos(incl) + y0
    
    dy_dx = np.gradient(y_vals, x_vals)
    
#     cond = x_vals <= xmax
    
#     x_vals = x_vals[cond]
#     y_vals = y_vals[cond]
#     dy_dx  = dy_dx[cond]

    for radius in range(5, int(rval_max), radial_steps):

        plt.clf()
        fig = plt.gcf()

        fig.set_figwidth(4)
        fig.set_figheight(4)

        axes_coords = [0, 0, 1, 1]

        ax = fig.add_axes(axes_coords)

        #im = ax.imshow(fits_file_obj.model.data, cmap='gray', alpha = 1, norm=colors.PowerNorm(0.3))
        im = ax.imshow(model_obj.data, cmap='gray', alpha = 1, norm=colors.PowerNorm(0.3))
        #fig.colorbar(im, orientation='vertical')
        ax.plot(x_vals, y_vals, color='springgreen', alpha=0.3,  linewidth=2)
        # Re-enable this comment to plot the arm opposite (assuming symmetry)
        #ax.plot(x_vals_flip, y_vals_flip, color='springgreen', alpha=0.3,  linewidth=2)

        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        #ax.set_ylim(-ymax/2, ymax/2)

        closest_idx = np.argmin(np.abs(rgrid - radius))
        closest_xval = x_vals[closest_idx]
        closest_yval = y_vals[closest_idx]
        closest_slope = dy_dx[closest_idx]
        
        closest_rval = rgrid[closest_idx]
        closest_theta = thetas[closest_idx]
        
        # # Plot circle
        circle1 = plt.Circle((x0, y0), closest_rval, fill = False, color = "yellow", linewidth = 1)
        ax.add_patch(circle1)

        xgrid1 = np.linspace(-xmax, closest_xval, num_pts//2)
        xgrid2 = np.linspace(closest_xval, xmax, num_pts//2)
        
        # Plot tangent to spiral
        tangent_line1 = closest_slope*(xgrid1 - closest_xval) + closest_yval
        tangent_line2 = closest_slope*(xgrid2 - closest_xval) + closest_yval
        
        ax.plot(xgrid1, tangent_line1, linestyle = "--", color = "orange")
        ax.plot(xgrid2, tangent_line2, linestyle = "--", color = "orange", label='_nolegend_')

        # # Plot tangent to circle
        # Thanks to https://www.quora.com/How-do-you-find-the-tangent-line-of-a-circle
        tangent_to_circle1 = (closest_xval**2 + closest_yval**2 + closest_rval**2 - 2*xgrid1*(closest_xval - x0) - x0**2 - y0**2)/(2*closest_yval - 2*y0)
        tangent_to_circle2 = (closest_xval**2 + closest_yval**2 + closest_rval**2 - 2*xgrid2*(closest_xval - x0) - x0**2 - y0**2)/(2*closest_yval - 2*y0)

        # Close...
        ax.plot(xgrid1, tangent_to_circle1, linestyle = "--", color = "royalblue")
        ax.plot(xgrid2, tangent_to_circle2, linestyle = "--", color = "royalblue", label='_nolegend_')

        plt.legend(["Spiral", f"Radius {radius}", "Tangent to Spiral", "Tangent to Circle"], loc = "upper right")

        # Validation
        circle_slope = np.gradient(tangent_to_circle1, xgrid1)[0]
        diff_angle   = np.degrees(np.arctan((closest_slope - circle_slope)/(1 + closest_slope*circle_slope)))
        #print(f"Difference in angle between spiral and circle: {180 - np.degrees(abs(angle_spiral - angle_circle)):.2f}")
        # Turn off to see difference in standard output
        #print(f"Difference in pitch angle measures (between tangents and differential) at radius {radius}: {abs(abs(diff_angle) - pitch_angles[closest_idx]):.2f}")
        
        pa_str = f"Derived PA: {pitch_angles[closest_idx]:.2f}"
        #print(f"PA at radius {closest_rval:.1f}:", f"{pitch_angles[closest_idx]:.2f}")
        ax.text(0.05, 0.95, pa_str, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', color = "white")#, bbox=props)
        
        pa_str = f"Geometric PA: {abs(diff_angle):.2f}"
        #print(f"PA at radius {closest_rval:.1f}:", f"{pitch_angles[closest_idx]:.2f}")
        ax.text(0.05, 0.90, pa_str, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', color = "white")#, bbox=props)

        gname_folder = pj(validation_dir, gname)
        if not exists(gname_folder):
            os.mkdir(gname_folder)

        filename = f'{pj(gname_folder, gname)}_{radius}.png'
        plt.savefig(filename)
        
    return filename

def generate_galfit_model_no_inclination(gname, feedme, out_dir, model_dir = "tmp_galfit_models"):
    
    #temp_dir    = pj(out_dir, "galfit_models")
    if not exists(model_dir):
        os.mkdir(model_dir)
        
    feedme_path = pj(model_dir, f"{gname}_model.in")
    
    feedme.header.output_image.value = pj(model_dir, f"{gname}_model_out.fits")
    feedme.header.optimize.value     = 1
    feedme.arms.inclination.value    = 0
    feedme.path_to_feedme            = feedme_path
    
    feedme.to_file()
    
    _ = sp(f"{run_galfit} {feedme_path}")
    
    return feedme.header.output_image.value

# Recommend turning use_inner_outer_rad on for pictures otherwise this is more informative
def calculate_pa(gpath, in_dir, out_dir, num_pts = 500, scatter_plot = True, validation_plot = False, use_inner_outer_rad = False, model_dir = "tmp_galfit_models"):
    
    gname = os.path.basename(gpath)
    try:
        fits_file = OutputFits(
            pj(gpath, f"{gname}_galfit_out.fits"), 
            #load_default = False,
            # Telling the module that these are the components we're looking for
            disk_for_arms = Sersic(3),
            arms          = Power(3),
            fourier       = Fourier(3),
            sky           = Sky(4),
            sersic_order  = ["bulge", "disk", "disk_for_arms"]
        )
        print(gname) #, 'ok')
    except (AttributeError, FileNotFoundError, Exception):
        print(f"Something went wrong opening galaxy {gname}! Continuing...")
        #continue
        return None
                                                        
    fit_region = fits_file.feedme.header.region_to_fit

    # Three(+1) component fit
    try:
        disk = fits_file.feedme.disk_for_arms
        
    # Two(+1) component fit
    except AttributeError:
        disk = fits_file.feedme.disk
        
    q          = disk.axis_ratio.value

    arms       = fits_file.feedme.arms
    fourier    = fits_file.feedme.fourier

    #amplitudes = [v[0] for k,v in fourier.param_values.items() if k != "skip"] #fourier.amplitudes
    #phis       = [v[1] for k,v in fourier.param_values.items() if k != "skip"] #fourier.phase_angles

    xmin = 0
    ymin = 0
    xmax = fit_region.x2 - fit_region.x1 #outer_rad*np.cos(arms.cumul_rot)
    ymax = fit_region.y2 - fit_region.y1 #outer_rad*np.sin(arms.cumul_rot)
    
    xgrid = np.linspace(xmin, xmax, num_pts)
    x0    = disk.position.x - fit_region.x1 #0 #0.5*xmax
    
    ygrid = np.linspace(ymin, ymax, num_pts)
    y0    = disk.position.y - fit_region.y1 #0 #0.5*ymax
    
    rvals_grid = np.sqrt((xgrid - x0)**2 + ((ygrid - y0)/q)**2)
    
    galaxy_info = pd.read_csv(pj(out_dir, gname, f"{gname}.csv"))
    crop_rad = float(galaxy_info[" cropRad"].iloc[0])
    scale_fact_std = 2*crop_rad/256

    try:
        arc_info = pd.read_csv(pj(out_dir, gname, f"{gname}_arcs.csv"))
        inner_rad = scale_fact_std*min(arc_info.loc[0, "r_start"], arc_info.loc[1, "r_start"])
        outer_rad = scale_fact_std*max(arc_info.loc[0, "r_end"], arc_info.loc[1, "r_end"])
    except:
        print(f"Something went wrong reading arc info from {gname}_arcs.csv")
        #continue
        return None
        
    if use_inner_outer_rad:
        #cond  = (rvals_grid > inner_rad) & (rvals_grid <= outer_rad)
        cond  = (rvals_grid > 0) & (rvals_grid <= 1.5*outer_rad)

    #     inner_idx = np.argmin(np.abs(np.abs(rvals_grid) - inner_rad))
    #     outer_idx = np.argmin(np.abs(np.abs(rvals_grid) - outer_rad)) + 1

    #     if outer_idx >= len(rvals_grid):
    #         outer_idx = len(rvals_grid) - 1


    #     old_inner = rvals_grid[:inner_idx]
    #     old_outer = rvals_grid[outer_idx:]

    #     # Basic dynamic grid, get higher resolution between inner and outer arm radius according to sparcfire
    #     high_res_xgrid = np.linspace(xgrid[inner_idx], xgrid[outer_idx - 1], num_pts)
    #     high_res_ygrid = np.linspace(ygrid[inner_idx], ygrid[outer_idx - 1], num_pts)
    #     high_res_rgrid = np.linspace(rvals_grid[inner_idx], rvals_grid[outer_idx - 1], num_pts)

    #     xgrid      = np.concatenate([old_inner, high_res_xgrid, old_outer])
    #     ygrid      = np.concatenate([old_inner, high_res_ygrid, old_outer])
    #     rvals_grid = np.concatenate([old_inner, high_res_rgrid, old_outer])

    #     num_pts    = len(rvals_grid)

    else:
        inner_rad = 0
        outer_rad = np.sqrt(xmax**2 + ymax**2)
        
        cond  = (rvals_grid > 0) & (rvals_grid <= outer_rad)

    if not np.any(cond):
        print("There was an issue constraining the radial grid via inner and outer rad.")
        return None

    xgrid = np.linspace(xgrid[cond][0], xgrid[cond][-1], num_pts)
    ygrid = np.linspace(ygrid[cond][0], ygrid[cond][-1], num_pts)
    
    rvals_grid = np.sqrt((xgrid - x0)**2 + ((ygrid - y0)/q)**2)

    # With Fourier Modes
    if fourier.parameters["skip"]:
        rgrid = rvals_grid
    else:
        rgrid = np.multiply(rvals_grid, 
                            sum_fourier_modes(
                                fourier.parameters,
                                #amplitudes = amplitudes,
                                #phis  = phis,
                                x0    = x0,
                                y0    = y0,
                                q     = q,
                                xgrid = xgrid,
                                ygrid = ygrid
                            )
                           )
    
    inner_idx = np.argmin(np.abs(np.abs(rgrid) - inner_rad))
    outer_idx = np.argmin(np.abs(np.abs(rgrid) - outer_rad))

    ones = np.ones(np.shape(rgrid))
    #with_fourier = np.array(list(map(pitch_angle, rgrid, ones*arms.cumul_rot, ones*arms.outer_rad, ones*arms.powerlaw, ones*arms.inner_rad)))
    pitch_angles, thetas = pitch_angle(
        rgrid, 
        ones*arms.cumul_rot.value, 
        ones*arms.outer_rad.value, 
        ones*arms.powerlaw_index.value, 
        ones*arms.inner_rad.value
    )
    
    #pitch_angles = pitch_angles[cond]
    #thetas       = thetas[cond]
    
    if scatter_plot:
        if not exists(pj(out_dir, "scatter_plots")):
            os.mkdir(pj(out_dir, "scatter_plots"))
        
        _ = plot_scatter(
                        pitch_angles, 
                        #rvals_grid[cond],
                        rgrid,
                        galaxy_info, 
                        arc_info, 
                        inner_idx = inner_idx, 
                        outer_idx = outer_idx, 
                        scatter_dir = pj(out_dir, "scatter_plots")
                        )
        
    if validation_plot:
        if not exists(pj(out_dir, "validation_plots")):
            os.mkdir(pj(out_dir, "validation_plots"))
            
        #thetas_rot = thetas - np.radians(arms.sky_position_angle) + 0.5*np.pi #- np.radians(disk.position_angle)
        thetas_rot = thetas + np.radians(arms.sky_position_angle.value) + 0.5*np.pi + np.radians(disk.position_angle.value)
        #rgrid_incl = rgrid*np.cos(np.radians(arms.inclination))
        
        model_output = generate_galfit_model_no_inclination(gname, fits_file.feedme, out_dir, model_dir = model_dir)
        try:
            model_only_fits_file = FitsFile(model_output)
        except (AttributeError, FileNotFoundError, Exception):
            print(f"Something went wrong opening galaxy model {gname}! Continuing with original output...")
            model_only_fits_file = fits_file.model
            #continue
            #return None
        
        _ = plot_validation(
                            pitch_angles, 
                            rgrid,
                            thetas_rot, 
                            galaxy_info, 
                            arc_info,
                            fits_file_obj = fits_file,
                            model_obj     = model_only_fits_file,
                            # inner_idx = inner_idx, 
                            # outer_idx = outer_idx, 
                            validation_dir = pj(out_dir, "validation_plots")
                            )

#     inner_idx = np.argmin(np.abs(rvals_grid - inner_rad))
#     outer_idx = np.argmin(np.abs(rvals_grid - outer_rad)) + 1

#     if outer_idx >= len(pitch_angles):
#         outer_idx = len(pitch_angles) - 1

    #avg_constrained_pa = np.mean(pitch_angles[inner_idx : outer_idx])
        
    #return avg_constrained_pa, rvals_grid[inner_idx], rvals_grid[outer_idx]
    return pitch_angles, thetas, rgrid #, inner_idx, outer_idx #rgrid[inner_idx : outer_idx]
      
def main(in_dir, out_dir, num_pts = 500, scatter_plot = True, validation_plot = False, model_dir = "tmp_galfit_models"):
    
    gpaths = glob(pj(out_dir, "123*"))
    gnames = [os.path.basename(i) for i in gpaths]
    
    pa_inner_outer = Parallel(n_jobs = -2)(
                           delayed(calculate_pa)( 
                               gpath,
                               in_dir,
                               out_dir,
                               num_pts             = num_pts,
                               scatter_plot        = scatter_plot,
                               validation_plot     = validation_plot,
                               use_inner_outer_rad = True,
                               model_dir           = model_dir
                                                )
                           for gpath in gpaths
                                            )
    
    pa_info_dict = dict(zip(gnames, pa_inner_outer))
    #print(pa_info_dict)
    
    avg_constrained_pas = [np.mean(tup[0]) if tup else None for tup in pa_inner_outer]
    
    pitch_angle_df = pd.DataFrame()
    pitch_angle_df = pitch_angle_df.from_dict(pa_info_dict, orient = "index", columns = ["pitch_angle", "theta", "fourier_rgrid"]) #, "inner_idx", "outer_idx"])
    
    pitch_angle_df.loc[:, "avg_constrained_pa"] = avg_constrained_pas
        
    return pitch_angle_df

if __name__ == "__main__":
    
    in_dir  = "sparcfire-in"
    tmp_dir = "sparcfire-tmp"
    out_dir = "sparcfire-out"
    
    scatter_plot    = True
    validation_plot = True
    cleanup = True
    
    if len(sys.argv) == 3:
        in_dir = sys.argv[1]
        out_dir = sys.argv[2]
        
    elif len(sys.argv) >= 4:
        in_dir = sys.argv[1]
        # Not used but this is the usual command line input
        tmp_dir = sys.argv[2]
        out_dir = sys.argv[3]
        
        try:
            true_tup = ("1", "true", "y")
            scatter_plot    = True if sys.argv[4].lower() in true_tup else False
            validation_plot = True if sys.argv[5].lower() in true_tup else False
        except IndexError:
            pass
        
    if not exists(in_dir) or not exists(out_dir):
        raise(f"Cannot find input/output directories {in_dir} {out_dir}.")
        
    model_dir = pj(out_dir, "tmp_galfit_models")
    
    pitch_angle_info = main(in_dir, out_dir, scatter_plot = scatter_plot, validation_plot = validation_plot, model_dir = model_dir)
    
    filename = pj(out_dir, "pitch-angle_info.pkl")
    pitch_angle_info.to_pickle(filename)
    #pd.dump(pa_rgrid_theta, open(filename, 'wb'))
    if cleanup:
        sp(f"rm -rf {model_dir}")
    