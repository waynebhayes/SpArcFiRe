import os
from os.path import join as pj
from os.path import exists

import pandas as pd
import shutil
from collections import namedtuple as nt

from Functions.helper_functions import sp
from Functions.helper_functions import in_notebook

from ResultsAnalysis.Plotting.create_all_plots import create_all_plots

from ResultsAnalysis.Helpers.load_petromags import load_petromags
#from ResultsAnalysis.Helpers.grab_header_parameters_for_flux import grab_header_parameters_for_flux
from ResultsAnalysis.Helpers.get_total_galaxies import get_total_galaxies
from ResultsAnalysis.Helpers.fprint import fprint
from ResultsAnalysis.Helpers.vprint import vprint

from ResultsAnalysis.residual_analysis import residual_analysis
from ResultsAnalysis.combine_multi_run_results import combine_multi_run_results
from ResultsAnalysis.create_quantiles import create_quantiles

def results_analysis(
    run_path, 
    *basenames, 
    **kwargs
):
    # Set some path variables and things
    run_path = run_path
    
    if in_notebook():
        run_path = run_path.replace("ics-home", "portmanm")

    in_dir  = kwargs.get("in_dir", pj(run_path, "sparcfire-in"))
    out_dir = kwargs.get("tmp_dir", pj(run_path, "sparcfire-out"))
    tmp_dir = kwargs.get("out_dir", pj(run_path, "sparcfire-tmp"))

    base_output_image_dir = kwargs.get("output_image_dir", pj(run_path, "for_paper_images"))     
    output_image_dirs = {bname : pj(base_output_image_dir, bname) for bname in basenames}
        
    _ = [os.makedirs(idir) for idir in output_image_dirs.values() if not exists(idir)]
    
    method          = kwargs.get("method", "nmr_x_1-p")
    nmr             = "norm_masked_residual"
    color_band      = kwargs.get("color_band", "r")
    
    global TOTAL_GALAXIES
    TOTAL_GALAXIES = get_total_galaxies(in_dir = in_dir, out_dir = out_dir)
    
    petromag_df = pd.DataFrame()
    # For speed
    if not kwargs.get("do_not_load_petromags"):
        print("Loading petromag info...")
        gzoo_file   = pj("/home", "portmanm", "kelly_stuff", "Kelly-29k.tsv")
        gzoo_file   = kwargs.get("gzoo_file", gzoo_file)

        petromag_df = load_petromags(gzoo_file, color_band)
    
    # FUNCTIONS OPTIONS
    incl_by_eye   = kwargs.get("incl_by_eye", False)
    by_eye_subset = kwargs.get("by_eye_subset", False)
    write         = kwargs.get("write", False)
    show          = kwargs.get("show", False)
    interactive   = kwargs.get("interactive", False)
    print_latex   = kwargs.get("print_latex", True)
    copy_png      = kwargs.get("copy_png", True)
    
    # Getting ready
    all_results  = {}
    plot_options = kwargs.get("plot_options", {bname : {} for bname in basenames})
    
    # LOOPING THROUGH NAMES GIVEN FOR ANALYSIS
    for basename in basenames:
        # RESIDUAL ANALYSIS
        fprint(f"PERFORMING RESIDUAL ANALYSIS FOR {basename}")
        analysis_results  = residual_analysis(
            in_dir              = in_dir, 
            out_dir             = out_dir, 
            basename            = basename,
            method              = method,
            incl_by_eye         = incl_by_eye,
            by_eye_subset       = by_eye_subset,
            color_band          = color_band,
            petromag_df         = petromag_df,
            pa_cutoff_val       = kwargs.get("pa_cutoff_val", 10),
            residual_cutoff_val = kwargs.get("residual_cutoff_val", 0.5),
            alen_cutoff_val     = kwargs.get("alen_cutoff_val", 0.007)
        )
        
        # Collating
        all_results[basename] = analysis_results
        
        if write or show:
            # OUTPUTTING PLOTS
            fprint("CREATING PLOTS")
            _ = create_all_plots(
                analysis_results, 
                method, 
                basename, 
                output_image_dirs[basename],
                incl_by_eye = incl_by_eye,
                write       = write,
                show        = show,
                interactive = interactive,
                pa_cutoff_val       = kwargs.get("pa_cutoff_val", 10),
                residual_cutoff_val = kwargs.get("residual_cutoff_val", 0.5),
                #alen_cutoff_val     = kwargs.get("alen_cutoff_val", 0.007)
                **plot_options[basename]
                #xaxis_range_mag_hist = [10, 20],
                #yaxis_range_mag_hist = [0, 0.15]
            )

        if print_latex or copy_png:
            fprint("QUANTILING IMAGES FROM RESULTS")
            
            if incl_by_eye:
                print("... by eye")
                quantile_df = analysis_results.by_eye_success_df
                galaxy_set_q = create_quantiles(
                    out_dir, 
                    #basename, 
                    quantile_df,
                    method,
                    **kwargs
                    # print_latex = print_latex, 
                    # copy_png = copy_png
                )
            else:
                quantile_df = analysis_results.success_df

            fprint("QUANTILING IMAGES FROM RESULTS")
            galaxy_set_q = create_quantiles(
                out_dir, 
                #basename, 
                quantile_df,
                method,
                **kwargs
                # print_latex = print_latex, 
                # copy_png = copy_png
            )
        
        # Unfortunately have to do this after and have the user generate the pngs from here
        # in order to rerun create_quantiles
        if kwargs.get("prep_for_quantile", False):
            fprint("JUST KIDDING, EXTRACTING QUANTILED MODELS TO BE CONVERTED TO PNG")
            
            to_untar = ' '.join([f"./{gname}_galfit_out.fits" for gname in galaxy_set_q])
            tar_file = f"{pj(out_dir, basename, basename)}_galfits.tar.gz"
            sp(f"tar -xzvf {tar_file} --occurrence {to_untar}")

            _ = [shutil.move(f"{gname}_galfit_out.fits", f"{pj(out_dir, basename, basename)}_galfits")
                 for gname in galaxy_set_q
                ]
            
            print(f"Please generate the pngs corresponding with the fits in the {pj(out_dir, basename, basename)}_galfits directory.")
            print("You may then proceed to run the 'create_quantiles' function again with copy_png set to True.")
    
    if len(basenames) > 1:
        fprint("COMBINING RESULTS FROM ALL RUNS FED IN")
        combined = combine_multi_run_results(
            method,
            *all_results.values(),
            df_names       = basenames,
            total_galaxies = TOTAL_GALAXIES,
            incl_by_eye    = incl_by_eye,
            by_eye_subset  = by_eye_subset
        )
        
        prefixes = list(set([i.split("_")[0] for i in basenames]))
        if len(prefixes) == 1:
            new_basename = f"{prefixes[0]}_combined"
        else:
            new_basename = "combined"
        
        all_results[new_basename]       = combined
        output_image_dirs[new_basename] = pj(base_output_image_dir, new_basename)
        if not exists(output_image_dirs[new_basename]):
            os.makedirs(output_image_dirs[new_basename])
            
        if write or show:
            try:
                # OUTPUTTING PLOTS
                fprint("CREATING PLOTS")
                _ = create_all_plots(
                    combined, 
                    method, 
                    new_basename, 
                    output_image_dirs[new_basename], 
                    incl_by_eye = incl_by_eye,
                    write       = write,
                    show        = show,
                    interactive = interactive,
                    pa_cutoff_val       = kwargs.get("pa_cutoff_val", 10),
                    residual_cutoff_val = kwargs.get("residual_cutoff_val", 0.007),
                    **plot_options[new_basename]
                    #xaxis_range_mag_hist = [10, 20],
                    #yaxis_range_mag_hist = [0, 0.15]
                )
            except KeyError as ke:
                print(f"Were plot options specified with the correct combined basename, {new_basename}?")
                print("Proceeding without plot options.")
                _ = create_all_plots(
                    combined, 
                    method, 
                    new_basename, 
                    output_image_dirs[new_basename], 
                    incl_by_eye = incl_by_eye,
                    write       = write,
                    show        = show,
                    interactive = interactive,
                    pa_cutoff_val       = kwargs.get("pa_cutoff_val", 10),
                    residual_cutoff_val = kwargs.get("residual_cutoff_val", 0.007),
                    #**plot_options[new_basename]
                    #xaxis_range_mag_hist = [10, 20],
                    #yaxis_range_mag_hist = [0, 0.15]
                )

        if print_latex or copy_png:
            fprint("QUANTILING IMAGES FROM RESULTS")
            if incl_by_eye:
                print("... by eye")
                # Do it twice to have a by eye sample and a regular sample
                quantile_df = combined.by_eye_success_df
                galaxy_set_q = create_quantiles(
                    out_dir, 
                    #basename, 
                    quantile_df,
                    method,
                    **kwargs
                    # print_latex = print_latex, 
                    # copy_png = copy_png
                )
                
            else:
                quantile_df = combined.success_df

            
            galaxy_set_q = create_quantiles(
                out_dir, 
                #basename, 
                quantile_df,
                method,
                **kwargs
                # print_latex = print_latex, 
                # copy_png = copy_png
            )
        
        # Unfortunately have to do this after and have the user generate the pngs from here
        # in order to rerun create_quantiles
        if kwargs.get("prep_for_quantile", False):
            fprint("JUST KIDDING, EXTRACTING QUANTILED MODELS TO BE CONVERTED TO PNG")
            
            to_untar = ' '.join([f"./{gname}_galfit_out.fits" for gname in galaxy_set_q])
            tar_file = f"{pj(out_dir, new_basename, new_basename)}_galfits.tar.gz"
            sp(f"tar -xzvf {tar_file} --occurrence {to_untar}")

            _ = [shutil.move(f"{gname}_galfit_out.fits", f"{pj(out_dir, new_basename, new_basename)}_galfits")
                 for gname in galaxy_set_q
                ]
            
            print(f"Please generate the pngs corresponding with the fits in the {pj(out_dir, basename, basename)}_galfits directory.")
            print("You may then proceed to run the 'create_quantiles' function again with copy_png set to True.")
        
    fprint("DONE!!!")
    
    output_format = """
    combined things only if applicable
    {
         basename : results namedtuple (fields below), 
         "[basename prefix]_combined" : combined results namedtuple,
    }
    (if the basename prefix isn't shared between the runs, the key is simply "combined")
    
    results namedtuple: 
        full_df, 
        success_df, 
        not_success_df, 
        by_eye_success_df, 
        by_eye_not_success_df
    
    combined results namedtuple: 
        bool_df,
        full_df,
        success_df, 
        by_eye_success_df, 
        
    purposefully chosen to mimic results above for analysis flexibility
    """
    
    return all_results
    