import pandas as pd
import numpy as np
from os.path import join as pj
from os.path import exists
from collections import namedtuple as nt
from copy import deepcopy

from ResultsAnalysis.Helpers.vprint import vprint
from Functions.helper_functions import find_files
from ResultsAnalysis.Helpers.grab_header_parameters_for_flux import grab_header_parameters_for_flux
from ResultsAnalysis.Helpers.get_total_galaxies import get_total_galaxies

ALL_RESULTS_NT = nt("ALL_RESULTS_NT", ["full_df", "success_df", "not_success_df", "by_eye_success_df", "by_eye_not_success_df"])

#=========================================================================
#=========================================================================

def load_residual_df(
    out_dir, 
    basename,
    **kwargs
):
    
    method              = kwargs.get("method", "nmr_x_1-p")
    verbose             = kwargs.get("verbose", True)
    residual_cutoff_val = kwargs.get("residual_cutoff_val", 0.007)
    
    pickle_filename = pj(out_dir, basename, sorted(find_files(pj(out_dir, basename), f'{basename}_output_results*.pkl', "f"))[-1])
    
    residual_df  = pd.read_pickle(pickle_filename)
    # temp_df = deepcopy(residual_df)
    # Setting residual columns
    #residual_df["KS_P"] = 1 - residual_df["KS_P"]
    if method == "nmr_x_1-p":
        result_of_method = (1 - residual_df["KS_P"])*residual_df["NMR"]
    elif method == "nmr_neg_log":
        result_of_method = residual_df["NMR"]/-np.log(residual_df["KS_P"] + 1e-10)
    elif method == "W_quality":
        result_of_method = residual_df["KS_P"]/residual_df["W_NMR"]
    else:
        raise Exception(f"Method given: {method} is not a valid method (yet).")
    
    residual_df[method] = result_of_method
    
    # Valid meaning NMR was successfully calculated
    #cols_to_drop = [col for col in residual_df.columns if col.endswith("_sky_2")]
    #valid_spiral_df = residual_df.drop(columns = cols_to_drop).dropna()

    # rename sky_2 to sky_3 for non-spirals to be inline with everything else
    # this would be for potential comparison down the line
    cols_to_merge = [col for col in residual_df.columns if col.endswith("_sky_3") or col.endswith("_sky_4")]
    #_ = [residual_df[col].fillna(residual_df[f"{col[:-1]}2"], inplace = True) for col in cols_to_merge]
    cols_to_drop  = [col for col in residual_df.columns if col.endswith("_sky_2") or col.endswith("_sky_3")]#  + ["KS_STAT"]
    residual_df.drop(columns = cols_to_drop, inplace = True)
    
    if verbose:
        print(f"{len(residual_df)} galaxy models generated.")
        residual_cutoff = residual_df[method] <= residual_cutoff_val
        print(f"{sum(residual_cutoff)} models pass score cutoff.")
    
    input_file = "./fake_file.fake_file"
    count      = 0
    while not exists(input_file) and count < 100:
        gname      = residual_df.index[count]
        input_file = pj(out_dir, gname, f"{gname}_galfit_out.fits")
        count += 1
        
    exptime, mag_zpt = grab_header_parameters_for_flux(input_file)
        
    # This will obviously have to change if multiple spiral components are introduced
    # but who's crazy enough to do that(???)
    spiral_component_num = int([i for i in residual_df.columns if i.startswith("inner_rad_power")][0][-1])
    
    # Formula from GALFIT readme...
    # m_tot = -2.5*log_10*(F_tot/t_exp) + mag_zpt
    # =>
    # F_tot = t_exp*10^[-(m_tot - mag_zpt)/2.5]    
    arm_flux   = exptime*10**(
        -(residual_df[f"magnitude_sersic_{spiral_component_num}"] - mag_zpt)/2.5
    )

    all_flux = deepcopy(arm_flux)
    for i in range(spiral_component_num - 1, 0, -1):
        flux = exptime*10**(
            -(residual_df[f"magnitude_sersic_{i}"] - mag_zpt)/2.5
        )
            
        all_flux += flux
        
        if i == 1:
            bulge_flux = flux

    residual_df["bulge_to_total_flux_ratio"] = bulge_flux/all_flux
    residual_df["arm_to_total_flux_ratio"]   = arm_flux/all_flux
    
    return residual_df.sort_values(by = method)

#=========================================================================
#=========================================================================

def load_galaxy_csv(out_dir, basename, pre_post):
    
    field = " pa_alenWtd_avg_domChiralityOnly"
    # {basename}_ uneccessary because different *galfit* runs 
    # should have same sparcfire output
    fname = pj(out_dir, basename, f"{basename}_{pre_post}_galfit_galaxy.csv")
    sparc_output_csv = pd.read_csv(fname, #pj(out_dir, f"pre_galfit_galaxy.csv"),
                                       index_col = "name",
                                       on_bad_lines = "warn",
                                       usecols   = ["name", field], # , " iptSz"],
                                       #na_values = "NaN",
                                       #dtype     = {field : float} #, " iptSz" : str}#, "name" : str}
                                      )#.loc[:, field]
    #sparc_output_csv.index.name = None
    sparc_output_csv[field] = sparc_output_csv[field].astype(float)
    sparc_output_csv.index  = sparc_output_csv.index.map(str)
    #sparc_output_csv[" iptSz"] = sparc_output_csv[" iptSz"].str.extract(r"([0-9]+)").astype(float)

    #sparc_output_csv["pre_sign"] = np.sign(sparc_output_csv[field])
    sparc_output_csv.rename(columns = {field : f"galaxy_{pre_post}_pa"}, inplace = True)
    
    return sparc_output_csv

#=========================================================================
#=========================================================================

def load_galaxy_arcs_csv(out_dir, basename, pre_post, **kwargs):
    
    field_pa   = kwargs.get("field_pa"  , "pitch_angle")
    field_alen = kwargs.get("field_alen", "arc_length")
    name_col   = kwargs.get("name_col"  , "gxyName")

    fname = pj(out_dir, basename, f"{basename}_{pre_post}_galfit_galaxy_arcs.csv")
    sparc_output_arcs_csv = pd.read_csv(fname, 
                                       index_col = name_col,
                                       usecols   = [name_col, field_pa, field_alen],
                                       dtype     = {field_pa : float, field_alen : float} #, name_col : str}
                                      )#.loc[:, field]
    #sparc_output_csv.index.name = None
    sparc_output_arcs_csv.index = sparc_output_arcs_csv.index.map(str)

    # Filtering for pure circles and near circles
    sparc_output_arcs_csv = sparc_output_arcs_csv[abs(sparc_output_arcs_csv[field_pa ]) > 1]

    #sparc_output_arcs_csv = pd.concat([sparc_output_arcs_csv, pre_sparc_output_csv], axis = 1)
    #sparc_output_arcs_csv["sign"] = np.sign(sparc_output_arcs_csv[field])

    # Keeps only arms which align with dom chirality only
    # sparc_output_arcs_csv["check"] = [
    #     row["sign"] + pre_sparc_output_csv.loc[i, "pre_sign"] 
    #     if i in pre_sparc_output_csv.index 
    #     else None 
    #     for i, row in sparc_output_arcs_csv.iterrows()
    # ]

    #sparc_output_arcs_csv = sparc_output_arcs_csv[abs(sparc_output_arcs_csv.loc[:, "check"]) == 2].drop(columns = ["sign", "check"])
    sparc_output_arcs_top3 = sparc_output_arcs_csv.groupby(name_col).head(3).reset_index()
    sparc_output_arcs_top3[f"{pre_post}_sign"] = np.sign(sparc_output_arcs_top3.pitch_angle)

    dom_sign = np.sign(sparc_output_arcs_top3.groupby(name_col).sum()[f"{pre_post}_sign"])
    sparc_output_arcs_top3 = sparc_output_arcs_top3.join(dom_sign, rsuffix = "_dom", on = name_col)

    cond = sparc_output_arcs_top3[f"{pre_post}_sign_dom"] == sparc_output_arcs_top3[f"{pre_post}_sign"]
    sparc_output_arcs_top2 = sparc_output_arcs_top3[cond].groupby(name_col).head(2).reset_index().drop(columns = [f"{pre_post}_sign_dom", "index"])

    #pre_sparc_output_top2.rename(columns = {field : "pre_pa"}, inplace = True)
    #pre_sparc_output_csv.dropna(inplace=True)
    return sparc_output_arcs_top2

#=========================================================================
#=========================================================================

def prepare_arcs_output(sparc_output_arcs_top2, pre_post, **kwargs):
    
    field_pa   = kwargs.get("field_pa"  , "pitch_angle")
    field_alen = kwargs.get("field_alen", "arc_length")
    name_col   = kwargs.get("name_col"  , "gxyName")
    
    single_arm = sparc_output_arcs_top2[~sparc_output_arcs_top2.duplicated(name_col, keep = False)]
    single_arm.loc[:, field_pa] = 0
    #single_arm.loc[:, "arc_length"]  = 0

    filled_in = pd.concat([sparc_output_arcs_top2, single_arm], ignore_index = True)
    str_fill = [f"{pre_post}_pa1", f"{pre_post}_pa2"] * (len(filled_in) // 2)
    filled_in["temp1"] = str_fill

    str_fill = [f"{pre_post}_alen1", f"{pre_post}_alen2"] * (len(filled_in) // 2)
    filled_in["temp2"] = str_fill

    #filled_in = filled_in.reset_index().drop(columns = ["index"])
    sp_out = filled_in.pivot_table(index = name_col, columns = ["temp1", "temp2"], values = [field_pa, field_alen])

    sp_out = sp_out.droplevel(0, axis = 1).droplevel(0, axis = 1)
    sp_out.columns = [f'{pre_post}_alen1', f'{pre_post}_alen2', f'{pre_post}_pa1', f'{pre_post}_pa2']
    
    return sp_out

#=========================================================================
#=========================================================================

def before_after_galfit_comparison(all_sparc_out, pre_sparc_output_csv, post_sparc_output_csv):
    
    before_after_galfit_df = all_sparc_out.copy(deep = True) #.dropna() #full_df.dropna(subset = ["post_pa"])
    #before_after_galfit_df = before_after_galfit_df[np.sign(before_after_galfit_df.loc[:, "pre_pa"]) != np.sign(before_after_galfit_df.loc[:, "post_pa"])]

    before_after_galfit_df["chiral_agreement"] = np.sign(before_after_galfit_df["pre_pa1"]) == np.sign(before_after_galfit_df["post_pa1"])

    before_after_galfit_df["pre_pa1"]  = abs(before_after_galfit_df["pre_pa1"])
    before_after_galfit_df["pre_pa2"]  = abs(before_after_galfit_df["pre_pa2"])
    before_after_galfit_df["post_pa1"] = abs(before_after_galfit_df["post_pa1"])
    before_after_galfit_df["post_pa2"] = abs(before_after_galfit_df["post_pa2"])
    
    #before_after_galfit_df.fillna(90, inplace = True)

    before_after_galfit_df["1-1"] = abs(before_after_galfit_df["pre_pa1"] - before_after_galfit_df["post_pa1"])
    before_after_galfit_df["2-2"] = abs(before_after_galfit_df["pre_pa2"] - before_after_galfit_df["post_pa2"])
    before_after_galfit_df["1-2"] = abs(before_after_galfit_df["pre_pa1"] - before_after_galfit_df["post_pa2"])
    before_after_galfit_df["2-1"] = abs(before_after_galfit_df["pre_pa2"] - before_after_galfit_df["post_pa1"])
    
    min_diff_pa_idx = before_after_galfit_df[["1-1", "2-2", "1-2", "2-1"]].idxmin(axis = 1)#.reset_index(drop = True)
    max_diff_pa_idx = before_after_galfit_df[["1-1", "2-2", "1-2", "2-1"]].idxmax(axis = 1)#.reset_index(drop = True)
    
    min_pre_str  = "pre_pa"  + min_diff_pa_idx.str[0]
    min_post_str = "post_pa" + min_diff_pa_idx.str[-1]
    #max_pre_str  = "pre_pa"  + min_diff_pa_idx.str[0]
    #max_post_str = "post_pa" + min_diff_pa_idx.str[-1]
    
    before_after_galfit_df["min_pre"], before_after_galfit_df["min_post"]  = zip(*[
        (before_after_galfit_df.loc[gname, col_name_pre], before_after_galfit_df.loc[gname, col_name_post])
        if 
            isinstance(col_name_pre, str) and isinstance(col_name_post, str)
        else 
            (None, None)
        for (gname, col_name_pre), (_, col_name_post) in zip(min_pre_str.items(), min_post_str.items())
        
    ])

    #before_after_galfit_df["pa_diff1"], before_after_galfit_df["pa_diff2"] = zip(*before_after_galfit_df["best_diffs"])
    #before_after_galfit_df["max_arm_pa_diff"]    = abs(before_after_galfit_df["max_post"] - before_after_galfit_df["max_pre"])
    before_after_galfit_df["min_arm_pa_diff"] = abs(before_after_galfit_df["min_post"] - before_after_galfit_df["min_pre"])
    before_after_galfit_df["pa_diff_galaxy"]  = abs(post_sparc_output_csv["galaxy_post_pa"] - pre_sparc_output_csv["galaxy_pre_pa"])
    
    before_after_galfit_df["min_pa_diff"]     = before_after_galfit_df[["min_arm_pa_diff", "pa_diff_galaxy"]].min(axis = 1)

    #before_after_galfit_df["alen_ratio"] = post_sparc_output_csv[" iptSz"]*before_after_galfit_df[["pre_alen1", "pre_alen2"]].min(axis = 1)/(pre_sparc_output_csv[" iptSz"]*before_after_galfit_df[["post_alen1", "post_alen2"]].max(axis = 1))
    before_after_galfit_df["alen_ratio"] = before_after_galfit_df[["post_alen1", "post_alen2"]].min(axis = 1)/before_after_galfit_df[["post_alen1", "post_alen2"]].max(axis = 1)

    cols_to_remove = ['1-1', '2-2', '1-2', '2-1'] # 'mean-1122', 'mean-1221', 'min_diff', 'best_diffs']
    before_after_galfit_df = before_after_galfit_df.drop(columns = cols_to_remove) #before_after_galfit_df.columns[9:-4])

    return before_after_galfit_df


def gather_everything(residual_df, before_after_galfit_df, method):
    full_df = residual_df.join(before_after_galfit_df)
    full_df = full_df[full_df.index.notnull()].sort_values(by = method)

    #full_df.dropna(subset = ["pa_diff1", "pa_diff2", "pa_diff_galaxy"], how = "all", inplace = True)
    #full_df.fillna(subset = ["pa_diff1", "pa_diff2", "pa_diff_galaxy"], how = "all", inplace = True)
    #full_df["min_pa_diff"] = full_df[["min_arm_pa_diff", "pa_diff_galaxy"]].min(axis = 1)
    
    return full_df

#=========================================================================
#=========================================================================

def determine_success(
    full_df, 
    **kwargs
):
    
    in_dir                = kwargs.get("in_dir", "sparcfire-in") 
    out_dir               = kwargs.get("out_dir","sparcfire-out")
    sparcfire_processed   = kwargs.get("sparcfire_processed", None)
    flip_chiral_agreement = kwargs.get("flip_chiral_agreement", False)
    residual_cutoff_val   = kwargs.get("residual_cutoff_val", 0.007)
    pa_cutoff_val         = kwargs.get("pa_cutoff_val", 10)
    alen_cutoff_val       = kwargs.get("alen_cutoff_val", 0.5)
    verbose               = kwargs.get("verbose", True)
    
    residual_cutoff = full_df["nmr_x_1-p"] <= residual_cutoff_val
    #pa_cutoff = (full_df["pa_diff1"] < 10) | (full_df["pa_diff2"] < 10)
    pa_cutoff   = full_df["min_pa_diff"] < pa_cutoff_val
    alen_cutoff = full_df["alen_ratio"] > alen_cutoff_val #[True]
    sign_cutoff = full_df["chiral_agreement"].astype(bool)
    #flux_ratio_cutoff = full_df["arm_to_total_flux_ratio"] > 0.05
    
    if flip_chiral_agreement:
        sign_cutoff = ~sign_cutoff

    success_df     = full_df[residual_cutoff & pa_cutoff & alen_cutoff & sign_cutoff].copy()
    not_success_df = full_df[~(residual_cutoff & pa_cutoff & alen_cutoff & sign_cutoff)].copy()
    
    if verbose:
        # print(f"{len(full_df)} processed by sparcfire")
        # print(f"{sum(residual_cutoff)} pass score cutoff")
        print(f"{sum(pa_cutoff)} pass pitch angle cutoff")
        print(f"{sum(alen_cutoff)} pass arm length ratio cutoff")
        print(f"{sum(sign_cutoff)} pass chiral agreement")
        #print(f"{sum(flux_ratio_cutoff)} pass arm flux ratio cutoff")
        
        print(f"{len(success_df)} or {100*len(success_df)/len(full_df):.2f}% ({len(success_df)}/{len(full_df)}) succeed by SpArcFiRe+Score")
        if sparcfire_processed is not None:
            sparcfire_processed = full_df.dropna(subset = ["min_arm_pa_diff", "pa_diff_galaxy"], how = "all")
        
        print(f"{TOTAL_GALAXIES - len(sparcfire_processed)}/{TOTAL_GALAXIES} models failed reprocessing by SpArcFiRe")
        
        #print(f"Total success less 24% false positive -- {len(success_df)*.76:.0f}")
        #print(f"Total success less 24% false positive + 24% false negative -- {len(not_success_df)*0.24+len(success_df)*.76:.0f}")
        #print(f"Estimated total success % -- {100*(len(not_success_df)*0.24+len(success_df)*.76)/len(full_df):.0f}%")
    
    # cutoffs = {
    #     "residual_cutoff" : residual_cutoff, 
    #     "pa_cutoff"       : pa_cutoff, 
    #     "alen_cutoff"     : alen_cutoff, 
    #     "sign_cutoff"     : sign_cutoff
    # }
    return success_df, not_success_df # , cutoffs

#=========================================================================
#=========================================================================

def extract_by_eye_data(
    out_dir, 
    basename, 
    residual_df, 
    full_df,
    **kwargs
):
    
    sparcfire_processed = kwargs.get("sparcfire_processed", None)
    subset              = kwargs.get("subset", None)
    verbose             = kwargs.get("verbose", True)
    
    with open(f"{pj(out_dir, basename, basename)}_by-eye_success.txt", "r") as f:
        raw_by_eye_success_galaxies = [i.split("_")[0].strip() for i in f.readlines()]

    with open(f"{pj(out_dir, basename, basename)}_by-eye_not_success.txt", "r") as f:
        raw_by_eye_not_success_galaxies = [i.split("_")[0].strip() for i in f.readlines()]
        
    by_eye_success_galaxies = [i for i in raw_by_eye_success_galaxies if i in full_df.index]
    by_eye_not_success_galaxies = [i for i in raw_by_eye_not_success_galaxies if i in full_df.index]
    if sparcfire_processed is not None:
        sparcfire_processed = full_df.dropna(subset = ["min_arm_pa_diff", "pa_diff_galaxy"], how = "all")
    
    if verbose:
        total = len(residual_df)
        if subset:
            total = subset
            print(f"Working on a subset of {total} galaxies")
            
        align = len(f"{len(by_eye_success_galaxies)}/{len(raw_by_eye_success_galaxies)}")
        print(f"Number of *total* by eye successful galaxies")
        print(f"{len(raw_by_eye_success_galaxies):<{align}} => {len(raw_by_eye_success_galaxies)/total*100:.2f}%")
        print(f"Number of by eye successful galaxies that SpArcFiRe *could* process")
        by_eye_processed = [i for i in sparcfire_processed.index if i in raw_by_eye_success_galaxies]
        print(f"{len(by_eye_processed)}/{len(raw_by_eye_success_galaxies)} => {len(by_eye_processed)/len(raw_by_eye_success_galaxies)*100:.2f}%")
        
        print()
        
        align = len(f"{len(by_eye_not_success_galaxies)}/{len(raw_by_eye_not_success_galaxies)}")
        print(f"Number of *total* by eye not successful galaxies")
        print(f"{len(raw_by_eye_not_success_galaxies):<{align}} => {len(raw_by_eye_not_success_galaxies)/total*100:.2f}%")
        
        print(f"Number of by eye not successful galaxies that SpArcFiRe *could* process")
        by_eye_processed = [i for i in sparcfire_processed.index if i in raw_by_eye_not_success_galaxies]
        print(f"{len(by_eye_processed)}/{len(raw_by_eye_not_success_galaxies)} => {len(by_eye_processed)/len(raw_by_eye_not_success_galaxies)*100:.2f}%")
    
    return by_eye_success_galaxies, by_eye_not_success_galaxies

#=========================================================================
#=========================================================================

def calculate_false_positive_negative(
    by_eye_success_galaxies, 
    by_eye_not_success_galaxies, 
    success_df, 
    not_success_df, 
    full_df,
    method  = "nmr_x_1-p",
    verbose = True
):
    
    false_positive = set(by_eye_not_success_galaxies).intersection(set(success_df.index))
    false_negative = set(by_eye_success_galaxies).intersection(set(not_success_df.index))

    by_eye_success_df     = full_df.loc[by_eye_success_galaxies].sort_values(by = method)
    by_eye_not_success_df = full_df.loc[by_eye_not_success_galaxies].sort_values(by = method)

    FP_rate = f"{len(false_positive)}/({len(false_positive)} + {len(by_eye_not_success_df)})"
    FN_rate = f"{len(false_negative)}/({len(false_negative)} + {len(by_eye_success_df)})"

    if verbose:
        print(f"False positive rate (by eye) -- {FP_rate} = {100*eval(FP_rate):.2f}%")
        print(f"False negative rate (by eye) -- {FN_rate} = {100*eval(FN_rate):.2f}%")

    #print(f"Total # of galaxies sorted by eye -- {len(raw_by_eye_success_galaxies) + len(raw_by_eye_not_success_galaxies)}")
    return by_eye_success_df, by_eye_not_success_df, FP_rate, FN_rate

#=========================================================================
#=========================================================================

def load_full_df_with_petromags(
    full_df,
    petromag_df,
    color_band
):
    
    petromag_col = f"petroMag_{color_band}"
        
    full_df   = full_df.join(petromag_df[petromag_col])#.dropna(subset = petromag_col)

    highest_mag_num = int(sorted([x for x in full_df.columns if x.startswith("magnitude_sersic")])[-1][-1])

    for i in range(1, highest_mag_num + 1):
        full_df[f"petromag_{color_band}_diff_{i}"] = full_df[f"magnitude_sersic_{i}"] - full_df[petromag_col]
        
    return full_df

#=========================================================================
#=========================================================================

def residual_analysis(
    **kwargs
):
    
    in_dir                = kwargs.get("in_dir", "sparcfire-in")
    out_dir               = kwargs.get("out_dir", "sparcfire-out")
    basename              = kwargs.get("basename", "") 
    method                = kwargs.get("method", "nmr_x_1-p")
    flip_chiral_agreement = kwargs.get("flip_chiral_agreement", False)
    pa_cutoff_val         = kwargs.get("pa_cutoff_val", 10)
    residual_cutoff_val   = kwargs.get("residual_cutoff_val", 0.5)
    alen_cutoff_val       = kwargs.get("alen_cutoff_val", 0.007)
    incl_by_eye           = kwargs.get("incl_by_eye", True)
    by_eye_subset         = kwargs.get("by_eye_subset", None)
    color_band            = kwargs.get("color_band", "r")
    petromag_df           = kwargs.get("petromag_df", pd.DataFrame())
    verbose               = kwargs.get("verbose", False)
    
    global TOTAL_GALAXIES
    TOTAL_GALAXIES = get_total_galaxies(in_dir = in_dir, out_dir = out_dir)
    
    vprint(verbose, "Load residual.")
    residual_df = load_residual_df(
        out_dir, 
        basename, 
        method = method, 
        residual_cutoff_val = residual_cutoff_val
    )
    
    # field_pa   = "pitch_angle"
    # field_alen = "arc_length"
    # name_col   = "gxyName"

    vprint(verbose, "Load pre galaxy csv.")
    pre_sparc_output_csv        = load_galaxy_csv(out_dir,      basename, pre_post = "pre")
        
    vprint(verbose, "Load pre galaxy arcs csv.")
    pre_sparc_output_arcs_top2  = load_galaxy_arcs_csv(out_dir, basename, pre_post = "pre")

    vprint(verbose, "Load post galaxy csv.")
    post_sparc_output_csv       = load_galaxy_csv(out_dir,      basename, pre_post = "post")
    
    vprint(verbose, "Load post galaxy arcs csv.")
    post_sparc_output_arcs_top2 = load_galaxy_arcs_csv(out_dir, basename, pre_post = "post")

#=========================================================================

    vprint(verbose, "Prep pre galaxy arcs df")
    pre_sp_out    = prepare_arcs_output(pre_sparc_output_arcs_top2,  pre_post = "pre")
    vprint(verbose, "Prep post galaxy arcs df")
    post_sp_out   = prepare_arcs_output(post_sparc_output_arcs_top2, pre_post = "post")

    vprint(verbose, "And combine")
    all_sparc_out = pd.concat([pre_sp_out, post_sp_out], axis = 1)
    
#=========================================================================

    vprint(verbose, "Compare SpArcFiRe analysis before and after")
    before_after_galfit_df     = before_after_galfit_comparison(
        all_sparc_out, 
        pre_sparc_output_csv, 
        post_sparc_output_csv
    )
    
    vprint(verbose, "Bring everything together")
    full_df             = gather_everything(residual_df, before_after_galfit_df, method)
    
    if not petromag_df.empty:
        vprint(verbose, "Load dataframe with petromag information")
        full_df         = load_full_df_with_petromags(full_df, petromag_df, color_band)
    
    sparcfire_processed = full_df.dropna(subset = ["min_arm_pa_diff", "pa_diff_galaxy"], how = "all")
    
    vprint(verbose, "Determine success")
    success_df, not_success_df = determine_success(
        full_df, 
        in_dir                 = in_dir, 
        out_dir                = out_dir, 
        flip_chiral_agreement  = flip_chiral_agreement,
        sparcfire_processed    = sparcfire_processed,
        pa_cutoff_val          = pa_cutoff_val, 
        residual_cutoff_val    = residual_cutoff_val,
        alen_cutoff_val        = alen_cutoff_val
    )
    
    full_df["success"] = full_df.index.isin(success_df.index)
    print()
    
#=========================================================================
    
    by_eye_success_df     = pd.DataFrame()
    by_eye_not_success_df = pd.DataFrame()
    
    if incl_by_eye:
        vprint(verbose, "Extract by-eye evaluation")
        by_eye_success_galaxies, by_eye_not_success_galaxies = extract_by_eye_data(
            out_dir, 
            basename, 
            residual_df, 
            full_df, 
            subset = by_eye_subset,
            sparcfire_processed = sparcfire_processed
        )
        print()

        # To resolve an occasional processing error...
        by_eye_success_limited     = list(set(by_eye_success_galaxies).intersection(full_df.index))
        by_eye_not_success_limited = list(set(by_eye_not_success_galaxies).intersection(full_df.index))

        vprint(verbose, "Calculate by-eye statistics")
        by_eye_success_df, by_eye_not_success_df, FP_rate, FN_rate = calculate_false_positive_negative(
            by_eye_success_limited, 
            by_eye_not_success_limited, 
            success_df, 
            not_success_df, 
            full_df,
            method = method
        )

        full_df["by_eye_success"] = full_df.index.isin(by_eye_success_df.index)
    
    results_nt = ALL_RESULTS_NT(full_df, success_df, not_success_df, by_eye_success_df, by_eye_not_success_df)
    for df in results_nt:
        df["runname"]  = basename
        
    return results_nt
