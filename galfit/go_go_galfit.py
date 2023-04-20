from galfit_objects import *
from sparc_to_galfit_feedme_gen import *
import shutil
import subprocess
# from os.path import join as pj
# from os.path import exists

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

def check_update(log_out, galfit_num, components, base_galfit_cmd, rerun = ""):
    
    #run_galfit, _, _ = check_programs()
    #base_galfit_cmd = f"{run_galfit} -imax {max_it}"
    
    # Put this here (for now) so it's easier to examine other components
    # of CompletedProcess from subprocess should I desire
    # i.e. I don't have to find it somewhere below
    log_out = log_out.stdout
    # Race condition issue
    #if exists(f"galfit.{galfit_num}"):
    
    # I don't like checking for the full line because there are embedded quotes
    # and one is a backtick I think...: Fit summary is now being saved into `fit.log'.
    # To be safe I just check the first part
    success = "Fit summary is now being saved"
    failure = "...now exiting to system..."
    
    # Constrain to last 30 lines to save search time
    last_out_lines = log_out.split("\n")[-30:]
    if any(line.startswith(success) for line in last_out_lines):
        final_model = log_out.split("Iteration")[-1]

        # bulge, disk, arms, fourier, sky
        components = update_components(final_model, 
                                       *components)
        
        galfit_num = f"{int(galfit_num) + 1:0>2}"
        
        if rerun:
            print("Rerunning GALFIT...")
            # TODO: This will run into issues with race conditions
            run_galfit_cmd = f"{base_galfit_cmd} galfit.{galfit_num}"
            galfit_num = f"{int(galfit_num) + 1:0>2}"
            out_text = sp(run_galfit_cmd)
            # Recursion... but forcing false to avoid rerunning over and over
            components, galfit_num = check_update(out_text, 
                                                  galfit_num, 
                                                  components, 
                                                  base_galfit_cmd, 
                                                  rerun = "")
            
    elif any(line.startswith(failure) for line in last_out_lines):
        print(f"Galfit failed this run!")
        # For debugging
        #print(log_out)
        #print("Expect output .# files to be less one albeit not necessarily missing this one.")
        
    else:
        print(f"Did not detect either '{success}' or '{failure}' in galfit output. Something must have gone terribly wrong! Printing output...")
        last_out_lines = '\n'.join(last_out_lines)
        print(f"{last_out_lines}")
        
    return components, galfit_num

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
    print(galaxy_names)
    if len(galaxy_names) == 1:
        gname = galaxy_names[0]
        
    print("Running feedme generator...")
    feedme_info = write_to_feedmes(top_dir = cwd, single_galaxy_name = gname)
    
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
        header  = feedme_info[gname]["header"]
        
        # These get *really* updated
        bulge   = feedme_info[gname]["bulge"]
        disk    = feedme_info[gname]["disk"]
        arms    = feedme_info[gname]["arms"]
        fourier = feedme_info[gname]["fourier"]
        sky     = feedme_info[gname]["sky"]
        
        components = (bulge, disk, arms, fourier, sky)
        
        galfit_num = "01"
        
        if num_steps == 1:
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_info[gname]['path']}"
            final_out_text = sp(run_galfit_cmd)
            
            # Dropping this here for rerun
            _, _ = check_update(final_out_text, 
                                galfit_num, 
                                components, 
                                base_galfit_cmd, 
                                rerun = rerun)
            
        elif num_steps >= 2:
            bulge_in = pj(out_dir, gname, f"{gname}_bulge.in")
            header.to_file(bulge_in, bulge, sky)
            
            run_galfit_cmd = f"{base_galfit_cmd} {bulge_in}"
            print("Bulge")
            out_text = sp(run_galfit_cmd)
            
            components, galfit_num = check_update(out_text, 
                                                  galfit_num, 
                                                  components, 
                                                  base_galfit_cmd, 
                                                  rerun = rerun)
            
            if num_steps == 3:
                disk_in = pj(out_dir, gname, f"{gname}_disk.in")
                header.to_file(disk_in, bulge, disk, sky)
                
                run_galfit_cmd = f"{base_galfit_cmd} {disk_in}"
                print("Bulge + Disk")
                out_text = sp(run_galfit_cmd)
                
                components, galfit_num = check_update(out_text, 
                                                      galfit_num, 
                                                      components, 
                                                      base_galfit_cmd, 
                                                      rerun = rerun)
                
            # Overwrite original... for now 
            header.to_file(f"{feedme_info[gname]['path']}", *components)
            run_galfit_cmd = f"{base_galfit_cmd} {feedme_info[gname]['path']}"
            print("Bulge + Disk + Arms")
            final_out_text = sp(run_galfit_cmd)
        
        tmp_png_path  = pj(tmp_png_dir, gname)
        tmp_fits_path = pj(tmp_fits_dir, f"{gname}_galfit_out.fits")
        
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
    