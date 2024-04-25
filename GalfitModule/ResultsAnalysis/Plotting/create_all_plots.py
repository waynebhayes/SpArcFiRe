import pandas as pd
import plotly.express as px

from ResultsAnalysis.Plotting.create_plot import create_plot

def create_all_plots(
    df_container,
    method, 
    basename, 
    output_image_dir,
    **kwargs
):
    
    #color_band   = kwargs.get("color_band", "r")
    incl_by_eye = kwargs.get("incl_by_eye", True)
    show        = kwargs.get("show", False)
    interactive = kwargs.get("interactive", False)
    write       = kwargs.get("write", False)

#=========================================================================
# FULL ECDF
#=========================================================================

    # Use a cutoff because there tends to be some extremely high values which skew the plot
    plot_df = df_container.full_df[df_container.full_df.loc[:, method] < kwargs.get("score_ecdf_cutoff", 0.015)].copy()
    
    _ = create_plot(
        x                = method,
        runname          = basename,
        plot_type        = "ecdf",
        df_or_dfs        = plot_df,
        output_image_dir = output_image_dir,
        xaxis_title      = method, #"KStest+NMR",
        marginal         = "histogram",
        # title       = f"1000 galaxies: ECDF for KStest+NMR on all models",
        # title_y     = 0.92
        cutoff_val       = kwargs.get("residual_cutoff_val", 0.007),
        show             = show,
        interactive      = interactive,
        write            = write
    )
    
#=========================================================================
# ECDF OF BY EYE SCORE
#=========================================================================
    
    if incl_by_eye:
        _ = create_plot(
            x                = method,
            runname          = f"{basename}_by-eye",
            plot_type        = "ecdf",
            df_or_dfs        = df_container.by_eye_success_df,
            output_image_dir = output_image_dir,
            xaxis_title      = method, #"KStest+NMR",
            marginal         = "histogram",
            #add_hline        = False,
            # title       = f"1000 galaxies: ECDF for KStest+NMR on by-eye successful model fits",
            # title_y     = 0.92
            show             = show,
            interactive      = interactive,
            write            = write
        )
    
#=========================================================================
# SERSIC INDEX HISTOGRAMS
#=========================================================================

    # x1   = "sersic_index_sersic_1"
    # x2   = "sersic_index_sersic_2"
    # x3   = "sersic_index_sersic_3"
    
    x    = "n"
    
    cols_to_use = sorted([x for x in df_container.full_df.columns if x.startswith("sersic_index_sersic")])        
    rename      = {col_name : f"sersic_{i}" for i, col_name in enumerate(cols_to_use)}
    
    plot_df_all          = df_container.full_df[cols_to_use].copy().rename(columns = rename)
    plot_df_success      = df_container.success_df[cols_to_use].copy().rename(columns = rename)
    
    to_plot  = {"all models" : plot_df_all, "success" : plot_df_success}
    runname   = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df_by_eye = df_container.by_eye_success_df[cols_to_use].copy().rename(columns = rename)
        to_plot["by-eye success"] = plot_df_by_eye
        runname   = f"{basename}_all-vs-success-vs-by-eye"
    
    _ = create_plot(
        x                       = x,
        runname                 = runname,
        plot_type               = "overlay_histogram",
        df_or_dfs               = to_plot,
        output_image_dir        = output_image_dir,
        histnorm                = "probability",
        #marginal                = "rug",
        #color                   = "component",
        #color_discrete_sequence = colors,
        nbins                   = kwargs.get("nbins", 40),
        #facet_col               = fcol,
        xaxis_range             = kwargs.get("xaxis_range_sersic_hist", None),
        yaxis_range             = kwargs.get("yaxis_range_sersic_hist", None),
        #reversed_order          = True,
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    
#=========================================================================
# MAGNITUDE HISTOGRAMS 
#=========================================================================

    # x1   = "magnitude_sersic_1"
    # x2   = "magnitude_sersic_2"
    # x3   = "magnitude_sersic_3"
        
    x    = "m"
        
    cols_to_use = sorted([x for x in df_container.full_df.columns if x.startswith("magnitude_sersic")])      
    rename      = {col_name : f"sersic_{i}" for i, col_name in enumerate(cols_to_use)}
    
    plot_df_all          = df_container.full_df[cols_to_use].copy().rename(columns = rename)
    plot_df_success      = df_container.success_df[cols_to_use].copy().rename(columns = rename)
    
    to_plot  = {"all models" : plot_df_all, "success" : plot_df_success}
    runname  = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df_by_eye            = df_container.by_eye_success_df[cols_to_use].copy().rename(columns = rename)
        to_plot["by-eye success"] = plot_df_by_eye
        runname                   = f"{basename}_all-vs-success-vs-by-eye"
    
    _ = create_plot(
        x                       = x,
        runname                 = runname,
        plot_type               = "overlay_histogram",
        df_or_dfs               = to_plot,
        output_image_dir        = output_image_dir,
        histnorm                = "probability",
        nbins                   = kwargs.get("mag_nbins"),
        hist_offset             = kwargs.get("mag_hist_offset", 3),
        xaxis_range             = kwargs.get("xaxis_range_mag_hist", None),
        yaxis_range             = kwargs.get("yaxis_range_mag_hist", None),
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    
#=========================================================================
# PETROMAG DIFFERENCE HISTOGRAM
#=========================================================================

    if "petromag" in " ".join(df_container.full_df.columns):
        x    = "m-petromag"
        fcol = "domain"

        cols_to_use = sorted([x for x in df_container.full_df.columns if x.startswith(f"petromag")])
        rename      = {col_name : f"sersic_{i}" for i, col_name in enumerate(cols_to_use)}
    
        plot_df_all          = df_container.full_df[cols_to_use].copy().rename(columns = rename)
        plot_df_success      = df_container.success_df[cols_to_use].copy().rename(columns = rename)

        to_plot  = {"all models" : plot_df_all, "success" : plot_df_success}
        runname  = f"{basename}_all-vs-success"

        if incl_by_eye:
            plot_df_by_eye            = df_container.by_eye_success_df[cols_to_use].copy().rename(columns = rename)
            to_plot["by-eye success"] = plot_df_by_eye
            runname                   = f"{basename}_all-vs-success-vs-by-eye"

        _ = create_plot(
            x                       = x,
            runname                 = runname,
            plot_type               = "overlay_histogram",
            df_or_dfs               = to_plot,
            output_image_dir        = output_image_dir,
            histnorm                = "probability",
            nbins                   = kwargs.get("mag_nbins"),
            hist_offset             = kwargs.get("mag_hist_offset", 3),
            xaxis_range             = kwargs.get("xaxis_range_petromag_hist", None),
            yaxis_range             = kwargs.get("yaxis_range_petromag_hist", None),
            # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
            # title_y     = 0.85
            show                    = show,
            interactive             = interactive,
            write                   = write
        )
        
#=========================================================================
# QUANTITIES SCATTER MATRIX
#=========================================================================

    #fcol  = "domain"
    #xmin, xmax = kwargs.get("xaxis_range_arm_flux_hist", (0, 5))
    to_compare  = ("magnitude", "effective", "sersic")
    #for name, df in zip(df_container._fields, df_container):
    #    if name not in ("full_df", "success_df", "by_eye_success_df"):
    #        continue        
    cols = df_container.full_df.columns
    #    name = name.replace("_df", "")
        
    all_index = [i for i in cols if i.startswith("sersic")]
    all_radii = [i for i in cols if i.startswith("effective")]
    all_mag   = [i for i in cols if i.startswith("magnitude")]
    
    plot_dfs = []
    # Limit sersic indices
    all_cutoffs  = [
        (df_container.full_df[n]   <= kwargs.get("index_cutoff",  4))
        &
        (df_container.full_df[re]  <= kwargs.get("radii_cutoff", 50))
        &
        (df_container.full_df[mag] <= kwargs.get("mag_cutoff",   20))
        &
        (df_container.full_df[mag] >= kwargs.get("mag_cutoff",   10))
        for n, re, mag in zip(all_index, all_radii, all_mag)
    ]

    df_cutoffs = all_cutoffs[0]
    for cutoff in all_cutoffs[1:]:
        # Boolean equivalent of and
        df_cutoffs *= cutoff

    plot_df = df_container.full_df[df_cutoffs].copy()
    plot_df["success_type"] = [
            2 if row["by_eye_success"]
        else
            1 if row["success"]
        else
            0
        for _, row in plot_df.iterrows()
    ]
    
#     # Avoid inplace just in case it resorts the original due to mutability(?)
    plot_df = plot_df.sort_values(by = ["success_type"])
    
    success_type_dict = {
        2 : "by eye success",
        1 : "by metric success",
        0 : "not successsful"
    }
    
    #We can do inplace here for sure
    plot_df["success_type"].replace(success_type_dict, inplace = True)
    
    color_discrete_map = {
        "by eye success" : px.colors.qualitative.Plotly[2],
        "success"        : px.colors.qualitative.Plotly[1],
        "not successful" : px.colors.qualitative.Plotly[0] 
    }
    
    color_order = ["not successful", "by metric success", "by eye success"]
    
    _ = create_plot(
        x                  = "n_re_m",
        runname            = basename,
        plot_type          = "scatter_matrix",
        df_or_dfs          = plot_df,
        # Magnitudes, effective radii, and sersic indices
        cols_to_compare    = [x for x in cols if x.split("_")[0] in to_compare],
        output_image_dir   = output_image_dir,
        color              = "success_type",
        #color_discrete_map = color_discrete_map,
        color_order        = color_order,
        height             = 1200,
        height_multiplier  = 1,
        #facet_col               = fcol,
        #xaxis_range             = kwargs.get("xaxis_range_scatter_matrix", None),
        #yaxis_range             = kwargs.get("yaxis_range_scatter_matrix", None),
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show               = show,
        interactive        = interactive,
        write              = write
    )
        
#=========================================================================
# SPIRAL ARM FLUX RATIO HISTOGRAM
#=========================================================================

    x     = "arm_to_total_flux_ratio"
    fcol  = "domain"
    #xmin, xmax = kwargs.get("xaxis_range_arm_flux_hist", (0, 5))

    plot_df       = df_container.full_df[x].copy().to_frame()
    plot_df[fcol] = "all models"

    plot_df1       = df_container.success_df[x].copy().to_frame()
    plot_df1[fcol] = "success"
    
    to_concat = [plot_df, plot_df1]
    runname   = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df2       = df_container.by_eye_success_df[x].copy().to_frame()
        plot_df2[fcol] = "by-eye success"
        to_concat.append(plot_df2)
        runname   = f"{basename}_by-eye-vs-success-vs-all"

    plot_df = pd.concat(to_concat, axis = 0)

    _ = create_plot(
        x                       = x,
        runname                 = runname,
        plot_type               = "histogram",
        df_or_dfs               = plot_df,
        output_image_dir        = output_image_dir,
        histnorm                = "probability",
        multi                   = False,
        facet_col               = fcol,
        nbins                   = kwargs.get("arm_flux_ratio_nbins"),
        #hist_offset             = kwargs.get("mag_hist_offset", 3),
        xaxis_range             = kwargs.get("xaxis_range_arm_flux_hist", None),
        yaxis_range             = kwargs.get("yaxis_range_arm_flux_hist", None),
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    
#=========================================================================
# B/T FLUX RATIO HISTOGRAM
#=========================================================================

    x     = "bulge_to_total_flux_ratio"
    fcol  = "domain"
    #xmin, xmax = kwargs.get("xaxis_range_arm_flux_hist", (0, 5))

    plot_df       = df_container.full_df[x].copy().to_frame()
    plot_df[fcol] = "all models"

    plot_df1       = df_container.success_df[x].copy().to_frame()
    plot_df1[fcol] = "success"
    
    to_concat = [plot_df, plot_df1]
    runname   = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df2       = df_container.by_eye_success_df[x].copy().to_frame()
        plot_df2[fcol] = "by-eye success"
        to_concat.append(plot_df2)
        runname   = f"{basename}_by-eye-vs-success-vs-all"

    plot_df = pd.concat(to_concat, axis = 0)

    _ = create_plot(
        x                       = x,
        runname                 = runname,
        plot_type               = "histogram",
        df_or_dfs               = plot_df,
        output_image_dir        = output_image_dir,
        histnorm                = "probability",
        multi                   = False,
        facet_col               = fcol,
        #nbins                   = 40, #kwargs.get("mag_nbins"),
        #hist_offset             = kwargs.get("mag_hist_offset", 3),
        xaxis_range             = kwargs.get("xaxis_range_bulge_flux_hist", None),
        yaxis_range             = kwargs.get("yaxis_range_bulge_flux_hist", None),
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    
#=========================================================================
# SPIRAL FLUX RATIO VS B/T FLUX RATIO SCATTER PLOTS
#=========================================================================

    x     = "bulge_to_total_flux_ratio"
    y     = "arm_to_total_flux_ratio"
    fcol  = "domain"
    #xmin, xmax = kwargs.get("xaxis_range_arm_flux_hist", (0, 5))

    plot_df       = df_container.full_df[x].copy().to_frame()
    plot_df[y]    = df_container.full_df[y].copy().to_frame()
    plot_df[fcol] = "all models"

    plot_df1       = df_container.success_df[x].copy().to_frame()
    plot_df1[y]    = df_container.success_df[y].copy().to_frame()
    plot_df1[fcol] = "success"
    
    to_concat = [plot_df, plot_df1]
    runname   = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df2       = df_container.by_eye_success_df[x].copy().to_frame()
        plot_df2[y]    = df_container.by_eye_success_df[y].copy().to_frame()
        plot_df2[fcol] = "by-eye success"
        to_concat.append(plot_df2)
        runname   = f"{basename}_by-eye-vs-success-vs-all"

    plot_df = pd.concat(to_concat, axis = 0)

    _ = create_plot(
        x                       = x,
        y                       = y,
        runname                 = runname,
        plot_type               = "scatter",
        df_or_dfs               = plot_df,
        output_image_dir        = output_image_dir,
        facet_col               = fcol,
        #nbins                   = 40, #kwargs.get("mag_nbins"),
        #hist_offset             = kwargs.get("mag_hist_offset", 3),
        xaxis_range             = kwargs.get("xaxis_range_flux_scatter", None),
        yaxis_range             = kwargs.get("yaxis_range_flux_scatter", None),
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    

#=========================================================================
# ALEN HISTOGRAM
#=========================================================================

    x    = "alen_ratio"
    fcol = "domain"
    
    plot_df       = df_container.full_df[x].copy().to_frame()
    plot_df[fcol] = "all models"

    plot_df1       = df_container.success_df[x].copy().to_frame()
    plot_df1[fcol] = "success"
    
    to_concat = [plot_df, plot_df1]
    runname   = f"{basename}_all-vs-success"
    
    if incl_by_eye:
        plot_df2       = df_container.by_eye_success_df[x].copy().to_frame()
        plot_df2[fcol] = "by-eye success"
        to_concat.append(plot_df2)
        runname   = f"{basename}_by-eye-vs-success-vs-all"

    plot_df = pd.concat(to_concat, axis = 0)
    
    _ = create_plot(
        x                       = x,
        runname                 = runname,
        plot_type               = "histogram",
        df_or_dfs               = plot_df,
        output_image_dir        = output_image_dir,
        histnorm                = "probability",
        multi                   = False,
        # color                   = "component",
        # color_discrete_sequence = colors,
        facet_col               = fcol,
        xaxis_range             = kwargs.get("xaxis_range_alen_hist"), # [10, 20],
        yaxis_range             = kwargs.get("yaxis_range_alen_hist"), # [0, 0.15],
        # title       = f"{runname} galaxies: distribution of magnitudes for by-eye successful models"
        # title_y     = 0.85
        show                    = show,
        interactive             = interactive,
        write                   = write
    )
    
#=========================================================================
# SCATTER OF PITCH ANGLE DIFFERENCES
#=========================================================================
    
    x     = "observation"
    y     = "model"
    color = "difference"

    plot_df        = df_container.full_df.copy()
    plot_df[x]     = df_container.full_df.loc[:, "min_pre"]
    plot_df[y]     = df_container.full_df.loc[:, "min_post"]
    plot_df[color] = df_container.full_df.loc[:, "min_arm_pa_diff"]

    _ = create_plot(
        x                = x,
        y                = y,
        runname          = f"{basename}_pa_diff",
        plot_type        = "scatter",
        df_or_dfs        = plot_df,
        output_image_dir = output_image_dir,
        color            = color,
        #xaxis_title      = "
        #color_continuous_midpoint = kwargs.get("pa_cutoff_val"),
        #range_color      = [0, 90],
        height_multiplier= 1,
        # title     = "Pitch angle difference reported by SpArcFiRe, model vs observation"
        # title_y   = 0.85
        show             = show,
        write            = write
    )
    
#=========================================================================
# SCATTER OF PITCH ANGLE DIFFERENCES FOR BY EYE SUCCESSFUL
#=========================================================================
    if incl_by_eye:
        x     = "observation"
        y     = "model"
        color = "difference"
        
        plot_df        = df_container.by_eye_success_df.copy()
        plot_df[x]     = df_container.by_eye_success_df.loc[:, "min_pre"]
        plot_df[y]     = df_container.by_eye_success_df.loc[:, "min_post"]
        plot_df[color] = df_container.by_eye_success_df.loc[:, "min_arm_pa_diff"]

        _ = create_plot(
            x                = x,
            y                = y,
            runname          = f"{basename}_by-eye_pa_diff",
            plot_type        = "scatter",
            df_or_dfs        = plot_df,
            output_image_dir = output_image_dir,
            color            = color,
            #color_continuous_midpoint = kwargs.get("pa_cutoff_val"),
            height_multiplier= 1,
            # title     = "Pitch angle difference reported by SpArcFiRe, model vs observation"
            # title_y   = 0.85
            show             = show,
            interactive      = interactive,
            write            = write
        )
    
#=========================================================================
# ECDF OF PITCH ANGLE DIFFERENCES
#=========================================================================

    # _ = create_plot(
    #     x                = "min_pa_diff",
    #     runname          = basename,
    #     plot_type        = "ecdf",
    #     df_or_dfs        = df_container.full_df,
    #     output_image_dir = output_image_dir,
    #     xaxis_title      = "Pitch Angle Difference (deg)",
    #     cutoff_val       = kwargs.get("pa_cutoff_val", 10),
    #     marginal         = "histogram",
    #     # title       = f"ECDF of pitch angle difference reported by SpArcFiRe, model vs observation"
    #     # title_y     = 0.85
    #     show             = show,
    #     interactive      = interactive,
    #     write            = write
    # )