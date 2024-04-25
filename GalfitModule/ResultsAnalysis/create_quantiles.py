import os
import shutil
from Functions.helper_functions import sp

def create_quantiles(
    out_dir, 
    df, 
    method,
    **kwargs
):
    print_latex = kwargs.get("print_latex", True)
    copy_png    = kwargs.get("copy_png", False)
    
    # Just in case
    df.sort_values(by = method, inplace = True)

    # Expect that if there exists more than one runname,
    # then we're working with combined data
    runnames = list(set(df.runname))
    if len(runnames) > 1:
        prefixes = list(set([i.split("_")[0] for i in runnames]))
        if len(prefixes) == 1:
            runname = f"{prefixes[0]}_combined"
        else:
            runname = "combined"
    else:
        runname = runnames[0]
        
    success_dir = pj(out_dir, runname, f'{runname}_galfit_png')
    print_latex_file = pj(out_dir, runname, f"{runname}_for_latex.txt")
    
    if kwargs.get("incl_by_eye", None):
        success_dir = pj(out_dir, runname, f'{runname}_by_eye_galfit_png')
        print_latex_file = pj(out_dir, runname, f"{runname}_by_eye_for_latex.txt")
    
    if not exists(success_dir):
        os.makedirs(success_dir)
    
    quantile           = ["0", "20", "40", "60", "80"]
    quantiled_galaxies = []
    
    print_latex_all = []
    if print_latex:
        
        if exists(print_latex_file):
            print("Deleting old latex output file...")
            os.remove(print_latex_file)
            
        print(f"Writing latex to file {print_latex_file}")
    
    for q in quantile:
        #vprint(print_latex, f"{q} &")
        print_latex_all.append(f"{q} &")
        
        if copy_png:
            quantile_dir = pj(success_dir, f"{runname}_all_quantile", f"quantile_{q}")
            if exists(quantile_dir):
                shutil.rmtree(quantile_dir)
            os.makedirs(quantile_dir)

        interp_df = df[method][df[method] >= df[method].quantile(0.01*float(q), interpolation='lower')]
        for count, (index, value) in enumerate(interp_df.items()):
            #if count < 5:
            #    continue
            if count == 8:
                break

            gname = index
            #print(q, i)
            #vprint(print_latex, f"{initial_str}{gname + '_combined.png'}{end_str}")
            #latex_rname = runname
            copy_rname  = df.loc[index, "runname"]
            
            # Use runname here for combined runs
            temp_str    = f"images/{runname}/{runname}_all_quantile/quantile_"
            initial_str = f"    \includegraphics[height=0.18\\textheight]{{{temp_str}{q}/"
            
            end_str = "} &"
            if count == 7 or count == len(interp_df) - 1:
                end_str = "} \\\\"
                
            print_latex_all.append(f"{initial_str}{gname + '_combined.png'}{end_str}")

            if copy_png:
                png_dir = pj(out_dir, copy_rname, f'{copy_rname}_galfit_png')
                shutil.copy(pj(png_dir, f"{gname}_combined.png"), quantile_dir)

            quantiled_galaxies.append(gname)
                
            #sp(f"cp {pj(out_dir, 'by_eye_success', gname + '_combined.png')} {pj(success_dir, 'all_quantile', 'quantile_' + q)}")
            
    if print_latex:           
        with open(print_latex_file, "w") as plf:
            plf.write("\n".join(print_latex_all))
            plf.write("\n")
            
    if copy_png:
        # Tar it all up!
        sp(f"tar -czvf {pj(out_dir, runname, runname)}_all_quantile.tar.gz -C {success_dir} {runname}_all_quantile")
        
    return quantiled_galaxies