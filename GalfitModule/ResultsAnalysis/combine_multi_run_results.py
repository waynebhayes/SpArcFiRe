import pandas as pd
from copy import deepcopy
from collections import namedtuple as nt

COMBINED_RESULTS_NT = nt("COMBINED_RESULTS_NT", ["bool_df", "full_df", "success_df", "by_eye_success_df"])

MINI_SEP    = "\n" + 40*"=" + "\n"

def combine_multi_run_results(
    method, 
    *args,
    **kwargs
):
    
    df_names       = kwargs.get("df_names", [])
    total_galaxies = kwargs.get("total_galaxies", 1)
    incl_by_eye    = kwargs.get("incl_by_eye", True)
    by_eye_subset  = kwargs.get("by_eye_subset", None)
    verbose        = kwargs.get("verbose", True)
    
    print(f"Joining {len(args)} attempts...")
    primary_full_df = deepcopy(args[0].full_df)
    
    num_dfs = len(args)
    #alt_full_df     = deepcopy(args[1].full_df)
    #alt_full_df.rename(columns = {method : f"1_{method}"}, inplace = True)

    all_full_dfs = [primary_full_df]
    all_methods  = [method]
    all_columns  = []
    
    for i, arg in enumerate(args[1:]):
        alt_method = f"{i}_{method}"
        all_methods.append(alt_method)
        
        all_full_dfs.append(arg.full_df.rename(columns = {method : alt_method}))
        all_columns.append(set(arg.full_df.columns))
        
    shared_columns = list(set(primary_full_df.columns).intersection(*all_columns)) + ["gname"]
    #empty_list = [None]*max([len(df) for df in all_full_dfs])
    #empty_df = pd.DataFrame({col : None for col in shared_columns}) #.set_index("gname")
    
    # BY RESIDUAL
    #temp_bool_df = pd.concat([primary_full_df[method], alt_full_df[f"1_{method}"]], axis = 1)
    combined_bool_df = pd.concat([df[method] for df, method in zip(all_full_dfs, all_methods)], axis = 1)

    #combined_bool_df.drop(index = list(set(primary_full_df.index).difference(set(alt_full_df.index))), inplace = True)
    #temp_bool_df["minima"] = temp_bool_df.idxmin(axis = 1)
    combined_bool_df["minima"] = combined_bool_df.idxmin(axis = 1)
    
    #og_minima  = temp_bool_df.minima == method
    #alt_minima = temp_bool_df.minima == f"1_{method}"
    
    #og_success = temp_bool_df.index.isin(args[0].by_eye_success_df.index)
    #alt_success = temp_bool_df.index.isin(args[1].by_eye_success_df.index)
    #print(sum((og_minima & og_success) | (alt_minima & alt_success)))
        
    # By everything
    eval_str = " | ".join([f"all_full_dfs[{i}].success" for i in range(num_dfs)])
    # success_n | success_m
    combined_bool_df["by_sparcfire_score_success"] = eval(eval_str) #primary_full_df.success | alt_full_df.success # | combined_bool_df.residual_minima_success
    
    minima_conditions  = [combined_bool_df.minima == method for method in all_methods]
    success_conditions = [df.success for df in all_full_dfs]
    all_conditions     = zip(minima_conditions, success_conditions)
    
    list_o_conditions = [cond_set[0] & cond_set[1] for cond_set in all_conditions]
    eval_str = " | ".join([f"list_o_conditions[{i}]" for i in range(num_dfs)])
    # minima -> success_minima
    combined_bool_df["by_sparcfire_and_best_score_success"] = eval(eval_str)
    
    combined_bool_df["best_fit"] = combined_bool_df[combined_bool_df.by_sparcfire_score_success].minima.replace(
        {alt_method : name for alt_method, name in zip(all_methods, df_names)}
    )
    
    # combined_bool_df["best_fit_and_success"] = combined_bool_df[combined_bool_df.by_sparcfire_score_success].minima.replace(
    #     {alt_method : name for alt_method, name in zip(all_methods, df_names)}
    # )
    
    if incl_by_eye:
        
        # Use by eye success df to account for by eye subsets (and to shorten the array) rather than info in full_df
        by_eye_success_conditions = [combined_bool_df.index.isin(df.by_eye_success_df.index) for df in args]
        # Flatten
        all_by_eye_success_gnames     = list(set([gname for df in args for gname in df.by_eye_success_df.index]))
        #all_by_eye_not_success_gnames = list(set([gname for df in args for gname in df.by_eye_not_success_df.index]))
        all_by_eye_not_success_gnames = list(set([gname for df in args for gname in df.by_eye_not_success_df.index]))
        #by_sparcfire_success_cond = [combined_bool_df.by_sparcfire_success == method for method in all_methods]

        combined_bool_df["residual_success_by_eye"] = combined_bool_df.index.isin(all_by_eye_success_gnames) & combined_bool_df.by_sparcfire_score_success 
        
        all_conditions = zip(minima_conditions, by_eye_success_conditions)
        list_o_conditions = [cond_set[0] & cond_set[1] for cond_set in all_conditions]
        eval_str = " | ".join([f"list_o_conditions[{i}]" for i in range(num_dfs)])

        # How well does choosing the smallest residual score across all runs work in picking a successful fit
        # when compared with the by eye analysis?
        # minima -> (success_minima & by eye)
        combined_bool_df["residual_minima_success_by_eye"] = eval(eval_str)
        # print(sum((minima_conditions[0] & by_eye_success_conditions[0]) | (minima_conditions[1] & by_eye_success_conditions[1])))
        # print(sum(list_o_conditions[0] | list_o_conditions[1]))
        
        by_sparcfire_success_by_eye = combined_bool_df.index.isin(all_by_eye_success_gnames) & combined_bool_df.by_sparcfire_score_success
            
        # As with and including residual minima, but now include the sparcfire scoring
        # (minima -> [success_minima & by eye]) | ([success_m | success_n] & by eye)
        combined_bool_df["by_minima_or_sparcfire_success_by_eye"]  = by_sparcfire_success_by_eye | combined_bool_df.residual_minima_success_by_eye
        # Comment out this one because it's filtering both individually by eye rather than doing (minima | score) & by eye
        #combined_bool_df["by_minima_and_sparcfire_success_by_eye"] = by_sparcfire_success_by_eye & combined_bool_df.residual_minima_success_by_eye
    
        # TODO: Show % in both
        # by eye success for all labeled by df
        best_fit_str_dict = {m : f"df_{i}" for i, m in enumerate(all_methods)}
        combined_bool_df["best_fit_by_eye"] = None

        for gname, row in combined_bool_df.iterrows():
            best_method = [
                (m, full_df.loc[gname, "by_eye_success"]) 
                for m, full_df in zip(all_methods, all_full_dfs)
                if gname in full_df.index and full_df.loc[gname, "by_eye_success"]
            ]

            if len(best_method) > 1:
                best_method = [(row.minima, None)]

            elif not best_method:
                best_method = [(None, None)]

            if not combined_bool_df.loc[gname, "best_fit_by_eye"]:
                combined_bool_df.loc[gname, "best_fit_by_eye"] = best_fit_str_dict.get(best_method[0][0], None)

        #eval_str = " | ".join([f"all_full_dfs[{i}].by_eye_success" for i, _ in enumerate(all_full_dfs)])
        #combined_bool_df["by_eye_success"] = eval(eval_str) #primary_full_df.by_eye_success | alt_full_df.by_eye_success
        combined_bool_df["by_eye_success"] = False | combined_bool_df.best_fit_by_eye.str.contains("df")
    
    if verbose:
        print(f"Total success by combining SpArcFiRe + score: {sum(combined_bool_df.by_sparcfire_score_success)}/{total_galaxies}")
        print(f"i.e. success_n | success_m | ...")
        print()
        print(f"Total success by combining SpArcFiRe + best score: {sum(combined_bool_df.by_sparcfire_and_best_score_success)}/{total_galaxies}")
        print(f"i.e. minima -> success_minima")
        print(MINI_SEP)
        if incl_by_eye:
            print("Checking against the by eye determination...")
            _total_galaxies = total_galaxies
            if df_names:
                if len(df_names) != num_dfs: 
                    print("Length of dataframe names supplied should be equal to the number of dataframes supplied.")
                    print("Leaving current convention in the dataframe (df_0, df_1, ..., df_n)")
                else:
                    combined_bool_df["best_fit_by_eye"]   = combined_bool_df.best_fit_by_eye.replace({f"df_{i}" : name for i, name in enumerate(df_names)})
            
            if by_eye_subset:
                _total_galaxies = by_eye_subset
                print("Using a subset of galaxies for the by eye determination...")
                
            print(f"Total success by eye: {sum(combined_bool_df.by_eye_success)}/{_total_galaxies}")
            total_by_eye = sum(combined_bool_df.by_eye_success)
            print()
            print(f"By eye captured by either score: {sum(combined_bool_df.residual_success_by_eye)}/{total_by_eye}")
            print(f"i.e. (success_m | success_n | ...) & by eye")
            print()
            print(f"By eye captured by best score: {sum(combined_bool_df.residual_minima_success_by_eye)}/{total_by_eye}")
            print(f"i.e. minima -> (success_minima & by eye)")
            print()
            print(f"By eye captured by SpArcFiRe or choosing best score between the two runs: {sum(combined_bool_df.by_minima_or_sparcfire_success_by_eye)}/{total_by_eye}")
            print(f"i.e. (minima -> [success_minima & by eye]) | ([success_m | success_n | ...] & by eye)")
            #print(f"By eye captured by SpArcFiRe and choosing best score between the two runs: {sum(combined_bool_df.by_minima_and_sparcfire_success_by_eye)}/{total_by_eye}")
            print(MINI_SEP)

            bss  = set(combined_bool_df[combined_bool_df.by_sparcfire_score_success].index)
            #bss  = set(combined_bool_df[bss.isin(all_by_eye_success_gnames)].index)
            TP   = all_by_eye_success_gnames
            #TP   = set(combined_bool_df[combined_bool_df["by_eye_success"]].index)
            
            #bsns  = ~combined_bool_df.by_sparcfire_score_success
            bsns = combined_bool_df[~combined_bool_df.by_sparcfire_score_success].index
            #bsns = bsns.index
            #bsns = set(combined_bool_df[bsns.isin(all_by_eye_not_success_gnames)].index)
            
            # Exclude the ones found in the success galaxies because some runs may find success where the others didn't
            TN =  set(all_by_eye_not_success_gnames).difference(set(all_by_eye_success_gnames))
            assert len(TP) + len(TN) == _total_galaxies, f"True positive and true negative don't add up to {_total_galaxies}!"
            #TN   = combined_bool_df[~combined_bool_df["by_eye_success"]].index
            #TN   = set(TN[TN.isin(all_by_eye_not_success_gnames)])

            FP   = bss.intersection(TN)
            FN   = bsns.intersection(TP)

            sparc_positive = bss.intersection(TP)
            sparc_negative = bsns.intersection(TN)
            fraction = len(sparc_positive)/sum(combined_bool_df.by_eye_success)
            
            combined_bool_by_eye_not_success = ~combined_bool_df.by_eye_success
            denom = combined_bool_by_eye_not_success[combined_bool_by_eye_not_success.index.isin(all_by_eye_not_success_gnames)]
            neg_fraction = len(sparc_negative)/sum(denom)
            # FPR = FP/(FP + TN)
            # FNR = FN/(FN + TP)
            # TODO WTF
            print(f"By eye success found by SpArcFiRe + score:  {len(sparc_positive)}/{sum(combined_bool_df.by_eye_success)} = {100*fraction:.2f}%")
            print(f"By eye not success found by SpArcFiRe + score:  {len(sparc_negative)}/{sum(denom)} = {100*neg_fraction:.2f}%")
            
            FP_rate = f"{len(FP)} / ({len(FP)} + {len(TN)})"
            FN_rate = f"{len(FN)} / ({len(FN)} + {len(TP)})"

            print()
            print(f"False positive rate (by eye) -- {FP_rate} = {100*eval(FP_rate):.2f}%")
            print(f"False negative rate (by eye) -- {FN_rate} = {100*eval(FN_rate):.2f}%")
            
            # TODO: GENERATE CONFUSION MATRIX
            #print()
            #print(f"Confusion matrix")
            #print()
            #print()
    
    
    _ = [full_df.rename(columns = {alt_method : method}, inplace = True) for alt_method, full_df in zip(all_methods[1:], all_full_dfs[1:])]
    
    combined_full_df = pd.concat([full_df for full_df in all_full_dfs])
    
    combined_success_df = pd.concat(
                full_df.loc[combined_bool_df[combined_bool_df.best_fit == name].index, :] 
                for name, full_df in zip(df_names, all_full_dfs)
            )
            
    combined_by_eye_success_df = None
    if incl_by_eye:
        # Get index, i.e. galaxy name from choosing the best fit then feed that into the full_dfs in all_full_dfs via loc
        # to grab the row
        combined_by_eye_success_df = pd.concat(
            full_df.loc[combined_bool_df[combined_bool_df.best_fit_by_eye == name].index, :] 
            for name, full_df in zip(df_names, all_full_dfs)
        )
    
    return COMBINED_RESULTS_NT(combined_bool_df, combined_full_df, combined_success_df, combined_by_eye_success_df)