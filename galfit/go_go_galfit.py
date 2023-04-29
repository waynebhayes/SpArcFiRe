from sparc_to_galfit_feedme_gen import *
import shutil
import subprocess
# from os.path import join as pj
# from os.path import exists

import sys
import os

_HOME_DIR = os.path.expanduser("~")
sys.path.append(f'{_HOME_DIR}/GalfitModule')
from Objects.Components import *
from Objects.Containers import *
from Functions.HelperFunctions import *

def sp(cmd_str, capture_output = True):
    # Because it is a pain in the butt to call subprocess with all those commands every time
    return subprocess.run(cmd_str, capture_output = capture_output, text = True, shell = True, executable="/bin/bash")

def check_programs():

    # This seems to work in Python directly so I'm leaving it as-is
    # Checking galfit
    run_galfit = shutil.which("galfit")
    #run_galfit = response.stdout.strip()

    # Checking fitspng
    run_fitspng   = shutil.which("fitspng")

    # Checking exact python3 call
    run_python = shutil.which("python3")

    return run_galfit, run_fitspng, run_python

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
        galfit_output = OutputContainer(sp(run_galfit_cmd), **galfit_output.to_dict())
       
    else:
        print("Galfit failed!")

    return galfit_output

def main(**kwargs):
    
    # Main Directories
    cwd       = kwargs.get("cwd", os.getcwd())
    in_dir    = kwargs.get("in_dir", pj(cwd, "sparcfire-in"))
    tmp_dir   = kwargs.get("tmp_dir", pj(cwd, "sparcfire-tmp"))
    out_dir   = kwargs.get("out_dir", pj(cwd, "sparcfire-out"))
    
    # Should some directory structure change in the future
    tmp_fits_dir = kwargs.get("tmp_fits_dir", pj(tmp_dir, "galfits"))
    tmp_png_dir  = kwargs.get("tmp_png_dir", pj(tmp_dir, "galfit_png"))
    out_png_dir  = kwargs.get("out_png_dir", pj(out_dir, "galfit_png"))
    
    # Of course the important things
    num_steps = int(kwargs.get("num_steps", 2))
    rerun     = kwargs.get("rerun", "")
    
    # If slurm, this will be a single filename
    galaxy_names = kwargs.get("galaxy_names", [])
    if isinstance(galaxy_names, str):
        galaxy_names = [galaxy_names]
    #assert isinstance(galaxy_names, list), "input_filenames must be a list, even if it's a single galaxy."
    
    # This PITA brought to you by running things in the notebook environment
    run_galfit, run_fitspng, run_python = check_programs()
    run_galfit  = kwargs.get("run_galfit", run_galfit)
    run_fitspng = kwargs.get("run_fitspng", run_fitspng)
    run_python  = kwargs.get("run_python", run_python)
    
    # We could chunk this up for slurm but since Wayne's script 
    # automatically distributes processes, we don't have to worry about it
    # Also this way we can take in one less variable (slurm)
    gname = ""
    slurm = False
    if len(galaxy_names) == 1:
        gname = galaxy_names[0]
        slurm = True
        
    print("Running feedme generator...")
    feedme_info = write_to_feedmes(top_dir = cwd,
                                   #single_galaxy_name = gname,
                                   in_dir  = in_dir,
                                   tmp_dir = tmp_dir,
                                   out_dir = out_dir
                                   )
    
    max_it = 150
    base_galfit_cmd = f"{run_galfit} -imax {max_it}"
    # TODO: Move output sigma file as well, figure out race condition
    # base_galfit_cmd = f"{run_galfit} -imax {max_it} -outsig"
    
    failed = []
    
    # Working on non-slurm for now
    print()
    print(f"Galfitting with {num_steps} component steps... (stdout is captured):")
    for gname in galaxy_names:
        # For debugging
        # if gname != "1237667783896924233":
        #     continue
        print(gname)
        
        # This doesn't change
        header  = feedme_info[gname].header
        feedme_path = feedme_info[gname].path_to_feedme
        # feedme_info[gname] is a FeedmeContainer object
        
        initial_components = feedme_info[gname].extract_components()
               
        # Note to self, OutputContainer is *just* for handling output
        # It does not retain the input state fed into Galfit
        # And the input dict is *before* output but will be updated *by* output
        #initial_galfit = OutputContainer(subprocess.CompletedProcess("",0), **components_dict)
        
        if num_steps == 1:
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_info[gname]['path']}"
            final_galfit_output = OutputContainer(sp(run_galfit_cmd), **feedme_info[gname].to_dict())
            
        elif num_steps >= 2:
            bulge_in = pj(out_dir, gname, f"{gname}_bulge.in")
            header.to_file(bulge_in, initial_components.bulge, initial_components.sky)
            
            run_galfit_cmd = f"{base_galfit_cmd} {bulge_in}"
            print("Bulge")
            galfit_output = OutputContainer(sp(run_galfit_cmd), **feedme_info[gname].to_dict())
            
            if rerun:
                galfit_output = rerun_galfit(galfit_output,
                                             base_galfit_cmd,
                                             initial_components.bulge, initial_components.sky
                                            )

            if num_steps == 3:
                disk_in = pj(out_dir, gname, f"{gname}_disk.in")
                
                # Note initial components disk and galfit output the rest
                # Those are updated!
                header.to_file(disk_in, galfit_output.bulge, initial_components.disk, galfit_output.sky)
                
                run_galfit_cmd = f"{base_galfit_cmd} {disk_in}"
                print("Bulge + Disk")
                galfit_output = OutputContainer(sp(run_galfit_cmd), **galfit_output.to_dict())
                
                if rerun:
                    galfit_output = rerun_galfit(galfit_output,
                                                 base_galfit_cmd,
                                                 galfit_output.bulge, initial_components.disk, galfit_output.sky
                                                )
    
            # Overwrite original... for now 
            # Also good thing dicts retain order, this frequently comes up
            galfit_output.to_file()
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_path}"
            print("Bulge + Disk + Arms")
            final_galfit_output = OutputContainer(sp(run_galfit_cmd), **galfit_output.to_dict())
            
        # Dropping this here for final rerun for all num_steps       
        if rerun:
            _ = rerun_galfit(final_galfit_output, 
                             base_galfit_cmd,
                             *galfit_output.to_list()
                            )

        tmp_png_path  = pj(tmp_png_dir, gname)
        tmp_fits_path = pj(tmp_fits_dir, f"{gname}_galfit_out.fits")
        
        # Alternatively, if not final_galfit_output.success and not galfit_output.success:
        if not exists(f"{tmp_fits_path}"):
            print(f"{gname} completely failed!!!")
            # This file avoids conflating ones which han'vet run and ones which have nothing to show for it
            sp(f"touch failed_{tmp_fits_path}", capture_output = False)
            # for debugging
            #print(final_out_text)
            print()
            failed.append(gname)
            continue
        
        fitspng_param = "0.25,1" #1,150"
        
        fitspng_cmd1 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}.png {tmp_fits_path}[1]"
        fitspng_cmd2 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_out.png {tmp_fits_path}[2]"
        fitspng_cmd3 = f"{run_fitspng} -fr \"{fitspng_param}\" -o \
                        {tmp_png_path}_residual.png {tmp_fits_path}[3]"
        
        _ = sp(fitspng_cmd1)
        _ = sp(fitspng_cmd2)
        _ = sp(fitspng_cmd3)
        
        # Combining the images using ImageMagick
        montage_cmd = f"montage {tmp_png_path}.png \
                                {tmp_png_path}_out.png \
                                {tmp_png_path}_residual.png \
                                -tile 3x1 -geometry \"175x175+2+0<\" \
                                {pj(out_png_dir, gname)}_combined.png"
        
        _ = sp(montage_cmd)
        
        # No point in doing this with slurm because race conditions
        if not slurm:
            for galfit_out in glob.glob(pj(cwd, "galfit.*")):
                ext_num = galfit_out.split(".")[-1]
                shutil.move(galfit_out, pj(out_dir, gname, f"{gname}_galfit.{ext_num}"))

        shutil.copy2(tmp_fits_path, pj(out_dir, gname, f"{gname}_galfit_out.fits"))
        
        print()
        
    return failed
        
if __name__ == "__main__":
    
    # cwd, in_dir, tmp_dir, out_dir, num_steps, rerun, *galaxy_names
    cmd_line_input = sys.argv[1:]
    
    kwargs_input = {}
    
    for kw in cmd_line_input:
        k, v = kw.split("=")
        kwargs_input[k] = v.strip()

    _ = main(**kwargs_input)
    