import shutil
import sys
import os
from os.path import join as pj
# from os.path import exists
from astropy.io import fits

_HOME_DIR = os.path.expanduser("~")    

try:
    _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
    _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
except KeyError:
    # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
    # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
    _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
sys.path.append(_MODULE_DIR)
    
from Classes.Components import *
from Classes.Containers import *
from Functions.helper_functions import *
from sparc_to_galfit_feedme_gen import *
from Utilities.parallel_residual_calc import fill_objects #, parallel_wrapper

import star_removal.no_log_remove_stars_with_sextractor as remove_stars_with_sextractor

# This is in helper_functions
# def check_programs():

#     # This seems to work in Python directly so I'm leaving it as-is
#     # Checking galfit
#     run_galfit = shutil.which("galfit")
#     #run_galfit = response.stdout.strip()

#     # Checking fitspng
#     run_fitspng   = shutil.which("fitspng")

#     # Checking exact python3 call
#     run_python = shutil.which("python3")

#     return run_galfit, run_fitspng, run_python

def rerun_galfit(galfit_output, base_galfit_cmd, *args):
    
    if galfit_output.success:   
        # don't use this anymore
        #galfit_num = f"{int(galfit_output.galfit_num) + 1:0>2}"
        
        print("Re-running GALFIT...")
        # String comparison until I implement subtraction ;)
        # Use this to generically determine which components were optimized
        updated_components = (comp for comp, icomp in zip(galfit_output.to_tuple(), args) if str(comp) != str(icomp))

        # run_galfit_cmd = f"{base_galfit_cmd} galfit.{galfit_num}"
        galfit_output.to_file(*updated_components)
        
        run_galfit_cmd = f"{base_galfit_cmd} {galfit_output.path_to_feedme}"

        #galfit_num = f"{int(galfit_output.galfit_num) + 1:0>2}"
        galfit_output = OutputContainer(sp(run_galfit_cmd, capture_output = True), **galfit_output.to_dict())
       
    else:
        print("Galfit failed!")

    return galfit_output

def touch_failed_fits(gname, tmp_fits_dir):
    
    prefix = "failed_"
    #filepath = os.path.dirname(tmp_fits_dir)
    fail_fileout = pj(tmp_fits_dir, f"{prefix}{gname}_galfit_out.fits")
    sp(f"touch {fail_fileout}", capture_output = False)
    
def main(**kwargs):
    
    # Main Directories
    cwd       = kwargs.get("cwd"    , os.getcwd())
    in_dir    = kwargs.get("in_dir" , pj(cwd, "sparcfire-in"))
    tmp_dir   = kwargs.get("tmp_dir", pj(cwd, "sparcfire-tmp"))
    out_dir   = kwargs.get("out_dir", pj(cwd, "sparcfire-out"))
    
    # Should some directory structure change in the future
    tmp_fits_dir  = kwargs.get("tmp_fits_dir" , pj(tmp_dir, "galfits"))
    tmp_png_dir   = kwargs.get("tmp_png_dir"  , pj(tmp_dir, "galfit_png"))
    tmp_masks_dir = kwargs.get("tmp_masks_dir", pj(tmp_dir, "galfit_masks"))
    tmp_psf_dir   = kwargs.get("tmp_psf_dir"  , pj(tmp_dir, "psf_files"))
    
    out_png_dir   = kwargs.get("out_png_dir"  , pj(out_dir, "galfit_png"))
    
    # Of course the important things
    num_steps = int(kwargs.get("num_steps", 2))
    rerun     = kwargs.get("rerun", False)
    parallel  = kwargs.get("parallel", 1)
    
    # For verbosity, default to capturing output
    # Keep both for clarity
    verbose        = kwargs.get("verbose", False)
    capture_output = kwargs.get("capture_output", True)
    
    # For generating starmasks per galaxy
    generate_starmasks = kwargs.get("generate_starmasks", True)
    
    # To indicate running from a local tmp dir
    run_from_tmp       = kwargs.get("run_from_tmp", False)
    
    # To aggressively free up space
    aggressive_clean   = kwargs.get("aggressive_clean", False)
    
    # Feeding in as comma separated galaxies
    # If single galaxy, this returns a list containing just the one galaxy
    galaxy_names = kwargs.get("galaxy_names", [])
    # petromags         = kwargs.get("petromags", [])
    # bulge_axis_ratios = kwargs.get("bulge_axis_ratios", [])
    if isinstance(galaxy_names, str):
        galaxy_names = galaxy_names.split(",")
        # petromags    = petromags.split(",")
        # bulge_axis_ratios    = bulge_axis_ratios.split(",")
    #     galaxy_names = [galaxy_names]
    #assert isinstance(galaxy_names, list), "input_filenames must be a list, even if it's a single galaxy."
    
    # This PITA brought to you by running things in the notebook environment
    run_galfit, run_fitspng, run_python = check_programs()
    run_galfit  = kwargs.get("run_galfit" , run_galfit)
    run_fitspng = kwargs.get("run_fitspng", run_fitspng)
    run_python  = kwargs.get("run_python" , run_python)
    
    # We could chunk this up for parallel but since Wayne's script 
    # automatically distributes processes, we don't have to worry about it
    # Deprecated
    # gname = ""
    # parallel = False
    # if len(galaxy_names) == 1:
    #     #gname = galaxy_names[0]
    #     parallel = True
    if run_from_tmp:
        for k,v in kwargs.items():
            if "_dir" in k and k != "out_dir":
                #_ = sp(f"rm -rf {v}")
            #os.makedirs(v, exist_ok = True)
                _ = sp(f"mkdir -p {v}")
    
    if generate_starmasks:
        print("Generating Starmasks")
        
        # I think we have to change path for source extractor
        star_removal_path = pj(_MODULE_DIR, "star_removal")
        os.chdir(star_removal_path)
        remove_stars_with_sextractor.main(in_dir, tmp_masks_dir, galaxy_names)
        os.chdir(cwd)
        
    print("Running feedme generator...")
    feedme_info = write_to_feedmes(top_dir = cwd,
                                   galaxy_names = galaxy_names,
                                   in_dir  = in_dir,
                                   tmp_dir = tmp_dir,
                                   out_dir = out_dir
                                   # petromags = petromags,
                                   # bulge_axis_ratios = bulge_axis_ratios
                                  )
                                   
    max_it = 200
    base_galfit_cmd = f"{run_galfit} -imax {max_it}"
    # TODO: Move output sigma file as well, figure out race condition
    # base_galfit_cmd = f"{run_galfit} -imax {max_it} -outsig"
    
    failed = []
    
    # Working on non-parallel for now
    print()
    print(f"Galfitting with {num_steps} component steps... (stdout is captured):")
    for gname in galaxy_names:
        # For debugging
        # if gname != "1237667783896924233":
        #     continue
        
        if gname not in feedme_info:
            print(f"Feedme not generated for {gname}")
            touch_failed_fits(gname, tmp_fits_dir)
            failed.append(gname)
            continue
            
        print(gname)
        
        disk_axis_ratio = 0.6
        with fits.open(pj(in_dir, f"{gname}.fits")) as gf:
            try:
                #SURVEY = SDSS-r  DR7
                #color = gf[0].header["SURVEY"].split()[0][-1]
                spirality_str = f"spirality"
                p_spirality = float(gf[0].header.get("spirality", 0))
                # p_spirality, disk_axis_ratio
                if p_spirality >= 0.9: 
                    disk_axis_ratio = 0.4
                elif p_spirality >= 0.75:
                    disk_axis_ratio = 0.5
                    
            except KeyError:
                print(f"There is no spirality keyword in the header for {gname}.")
        
        # This doesn't change
        header  = feedme_info[gname].header
        feedme_path = feedme_info[gname].path_to_feedme
        # feedme_info[gname] is a FeedmeContainer object
        
        initial_components = feedme_info[gname].extract_components()
               
        # Note to self, OutputContainer is *just* for handling output
        # It does not retain the input state fed into Galfit
        # And the input dict is *before* output but will be updated *by* output
        #initial_galfit = OutputContainer(subprocess.CompletedProcess("",0), **components_dict)
        
        # Any galfit runs that feed into an output container *must* have capture_output = True or 
        # they won't update the new parameters
        if num_steps == 1:
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_path}"
            final_galfit_output = OutputContainer(sp(run_galfit_cmd), **feedme_info[gname].to_dict())
            
        elif num_steps >= 2:
            # This ends up being more like a disk fit
            disk_in = pj(out_dir, gname, f"{gname}_disk.in")
            header.to_file(disk_in, initial_components.disk, initial_components.sky)
            
            run_galfit_cmd = f"{base_galfit_cmd} {disk_in}"
            print("Disk")
            galfit_output = OutputContainer(sp(run_galfit_cmd), sersic_order = ["disk"], **feedme_info[gname].to_dict())
            
            # Assume rerun is to refine final fit
            # This will also be better for prototyping multiband fits
            # if rerun:
            #     galfit_output = rerun_galfit(galfit_output,
            #                                  base_galfit_cmd,
            #                                  # Pass in initial components here to generically
            #                                  # determine which were optimized on
            #                                  initial_components.disk, initial_components.sky
            #                                 )

            if num_steps == 3:
                # For fitting *just* the bulge + a few pixels
                # Trying Just disk, just bulge, then all together with arms
#                 bulge_header = deepcopy(feedme_info[gname].header)
#                 xcenter, ycenter = initial_components.bulge.position
#                 bulge_rad = 2*initial_components.bulge.effective_radius
#                 bulge_header.region_to_fit = (
#                                               round(xcenter - bulge_rad), 
#                                               round(xcenter + bulge_rad), 
#                                               round(ycenter - bulge_rad), 
#                                               round(ycenter + bulge_rad)
#                                              )
#                 bulge_header.param_fix["region_to_fit"] = f"{round(ycenter - bulge_rad)} {round(ycenter + bulge_rad)}"
#                 bulge_header.update_param_values()

#                 bulge_in = pj(out_dir, gname, f"{gname}_bulge.in")
#                 bulge_header.to_file(bulge_in, initial_components.bulge, galfit_output.sky)

                # Fit disk first, now to fit the bulge and refine the disk
                bulge_disk_in = pj(out_dir, gname, f"{gname}_bulge+disk.in")
                
                # Note initial components disk and galfit output the rest
                # Those are updated!
                header.to_file(bulge_disk_in, initial_components.bulge, galfit_output.disk, galfit_output.sky)
                
                run_galfit_cmd = f"{base_galfit_cmd} {bulge_in} {bulge_disk_in}"
                print("Bulge + Disk")
                galfit_output = OutputContainer(sp(run_galfit_cmd), **galfit_output.to_dict())
                
#                 if rerun:
#                     galfit_output = rerun_galfit(galfit_output,
#                                                  base_galfit_cmd,
#                                                  initial_components.bulge, galfit_output.disk, galfit_output.sky
#                                                 )
    
            # Overwrite original... for now 
            # Also good thing dicts retain order, this frequently comes up
            if galfit_output.arms.param_values.get("skip", 0):
                # By default includes the header
                print("Skipping Arms")
                galfit_output.to_file(galfit_output.bulge, galfit_output.disk, galfit_output.sky)
            else:
                galfit_output.disk.axis_ratio = disk_axis_ratio
                galfit_output.disk.update_param_values()
                #galfit_output.arms.param_fix["outer_rad"] = 1
                galfit_output.to_file()
                
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_path}"
            print("Bulge + Disk + Arms (if applicable)")
            final_galfit_output = OutputContainer(sp(run_galfit_cmd), **galfit_output.to_dict(), store_text = True)
            
        # Dropping this here for final rerun for all num_steps
        if rerun:
            _ = rerun_galfit(final_galfit_output, 
                             base_galfit_cmd,
                             *galfit_output.to_list()
                            )
        
        if verbose:
            print(str(final_galfit_output))

        tmp_png_path  = pj(tmp_png_dir, gname)
        tmp_fits_path_gname = pj(tmp_fits_dir, f"{gname}_galfit_out.fits")
        
        # Alternatively, if not final_galfit_output.success and not galfit_output.success:
        if not exists(f"{tmp_fits_path_gname}"):
            print(f"{gname} completely failed!!!")
            # This file avoids conflating ones which haven't run and ones which have nothing to show for it
            touch_failed_fits(gname, tmp_fits_dir)
            # for debugging
            #print(final_out_text)
            print()
            #failed.append(gname)
            continue
        
        fitspng_param = "0.25,1" #1,150"
        
        fitspng_cmd1 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}.png {tmp_fits_path_gname}[1]"
        fitspng_cmd2 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_out.png {tmp_fits_path_gname}[2]"
        fitspng_cmd3 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_residual.png {tmp_fits_path_gname}[3]"
        
        fitspng_out = sp(fitspng_cmd1)
        
        # Because this doesn't work on some of the clusters
        if "error" not in fitspng_out.stderr:
            
            _ = sp(fitspng_cmd2, capture_output = capture_output)
            _ = sp(fitspng_cmd3, capture_output = capture_output)

            # Combining the images using ImageMagick
            montage_cmd = f"montage {tmp_png_path}.png \
                                    {tmp_png_path}_out.png \
                                    {tmp_png_path}_residual.png \
                                    -tile 3x1 -geometry \"175x175+2+0<\" \
                                    {pj(out_png_dir, gname)}_combined.png"

            _ = sp(montage_cmd, capture_output = capture_output)
            # Drop a copy in the sparfire-out folder of each galaxy for ease of navigating/viewing
            shutil.copy2(f"{pj(out_png_dir, gname)}_combined.png", pj(out_dir, gname))
                
        else:
            print("Skipping fitspng conversion... there is likely a library (libcfitsio) issue.")

        # No point in doing this in parallel because race conditions
        if not parallel:
            for galfit_out in glob.glob(pj(cwd, "galfit.*")):
                ext_num = galfit_out.split(".")[-1]
                shutil.move(galfit_out, pj(out_dir, gname, f"{gname}_galfit.{ext_num}"))
                
        # Residual calculation, now done all at the same time! And added to the FITS header
        _, gname_nmr, gname_kstest = fill_objects(gname, 1, tmp_fits_dir, tmp_masks_dir)
        with fits.open(tmp_fits_path_gname, mode='update', output_verify='ignore') as hdul:
            hdul[2].header["NMR"] = (gname_nmr, "Norm of the masked residual")
            if gname_kstest:
                hdul[2].header["ks_p"] = (round(gname_kstest.pvalue, 4), "p value of kstest vs noise")
                hdul[2].header["ks_stat"] = (round(gname_kstest.statistic, 4), "statistic value of kstest vs noise")
            else:
                hdul[2].header["ks_p"] = (None, "p value of kstest vs noise")
                hdul[2].header["ks_stat"] = (None, "statistic value of kstest vs noise")
            # Flush done automatically in update mode
            #hdul.flush()

        shutil.copy2(tmp_fits_path_gname, pj(out_dir, gname, f"{gname}_galfit_out.fits"))
        
        # TODO(?): Output final GALFIT parameters/components to file via FITS header and
        # FitsHandler routines. This may somewhat useful for future comparison but
        # theoretically any comparison will likely take place in python and can
        # therefore use the FitsHandler routines at run time.
        
        print()
        
    if aggressive_clean:
        print("Aggressively cleaning... in, temp galfit output, and psf directories.")
        to_del_gal_in    = (pj(in_dir, f"{gname}.fits") for gname in galaxy_names
                            if exists(pj(in_dir, f"{gname}.fits"))
                           )

        to_del_tmp_fits  = (pj(tmp_fits_dir, f"{gname}_galfit_out.fits") for gname in galaxy_names
                            if exists(pj(tmp_fits_dir, f"{gname}_galfit_out.fits"))
                           )
        
        to_del_masks     = (pj(tmp_masks_dir, f"{gname}_star-rm.fits") for gname in galaxy_names
                            if exists(pj(tmp_masks_dir, f"{gname}_star-rm.fits"))
                           )
        
        to_del_psf_files = (pj(tmp_psf_dir, f"{gname}_psf.fits") for gname in galaxy_names
                            if exists(pj(tmp_psf_dir, f"{gname}_psf.fits"))
                           )
        
        #print(f"rm -rf {' '.join(to_del_in)} {' '.join(to_del_tmp_fits)} {' '.join(to_del_psf_files)}")
        _ = sp(f"rm -f {' '.join(to_del_gal_in)} {' '.join(to_del_tmp_fits)} {' '.join(to_del_masks)} {' '.join(to_del_psf_files)}",
               capture_output = capture_output
              )
    
    return failed
        
if __name__ == "__main__":
    
    # cwd, in_dir, tmp_dir, out_dir, num_steps, rerun, verbose, *galaxy_names
    cmd_line_input = sys.argv[1:]
    
    kwargs_input = {}
    
    for kw in cmd_line_input:
        k, v = kw.split("=")
        v = v.strip()
        if v.lower() == "false":
            v = False
        elif v.lower() == "true":
            v = True
            
        kwargs_input[k] = v

    _ = main(**kwargs_input)
    #write_failures(os.getcwd(), failures)
    
