import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from copy import deepcopy
from math import ceil

#import matplotlib as plt

#=========================================================================
#=========================================================================

def create_ecdf(
    x,
    df, 
    dict_o_kwargs
):

    fig = px.ecdf(
        df,
        x        = x,
        markers  = True, 
        lines    = False, 
        marginal = dict_o_kwargs.get("marginal"),
        ecdfnorm = None,
     ) 

    cutoff_val = dict_o_kwargs.get("cutoff_val")
    if dict_o_kwargs.get("add_vline"):
        fig.add_vline(x = cutoff_val, 
                      row = 1,
                      line_color = "cyan",
                      annotation_text= f"{cutoff_val}", 
                      annotation_position="bottom")

    if dict_o_kwargs.get("add_hline"):
        yval = sum(df.loc[:, x] < cutoff_val)
        fig.add_hline(y = yval, 
                      row = 1,
                      col = 1,
                      line_color = "magenta",
                      annotation_text=f"{yval}",
                      annotation_position="bottom left"
                     )
        
    return fig

#=========================================================================
#=========================================================================

def create_scatter(
    x, 
    y, 
    df, 
    dict_o_kwargs
):
    
    color_continuous_scale    = "Agsunset"
    #if dict_o_kwargs.get("color_continuous_midpoint"):
    #    color_continuous_scale = "Portland"
    
    fig = px.scatter(
        df, 
        x = x, 
        y = y,
        facet_col = dict_o_kwargs.get("facet_col"),
        facet_row = dict_o_kwargs.get("facet_row"),
        color                  = dict_o_kwargs.get("color"),
        color_continuous_scale = color_continuous_scale,
        range_color            = dict_o_kwargs.get("range_color"),
        #color_continuous_midpoint = dict_o_kwargs.get("color_continuous_midpoint"),
    )
    
    if dict_o_kwargs.get("facet_col") or dict_o_kwargs.get("facet_row"):
        fig.for_each_annotation(lambda a: a.update(text = a.text.split("=")[-1]))
        
    return fig

#=========================================================================
#=========================================================================

# import seaborn as sns
# sns.set_theme(style="ticks")
def create_scatter_matrix(
    df,
    cols_to_compare,
    dict_o_kwargs,
    **kwargs
):
    
    #color_continuous_scale    = "Agsunset"
    #if dict_o_kwargs.get("color_continuous_midpoint"):
    #    color_continuous_scale = "Portland"
    
    rename_dict = {
        "magnitude" : "m",
        "effective" : "r_e",
        "sersic"    : "n",
    }
    
    labels = {}
    for col in cols_to_compare:
        num = col[-1]
        labels[col] = rename_dict.get(col.split('_')[0]) + num
    
#     fig = sns.pairplot(
#         df.replace(labels), 
#         hue       = dict_o_kwargs.get("color"),
#         hue_order = dict_o_kwargs.get("color_order"),
#         palette   = dict_o_kwargs.get("color_discrete_map"),
#         vars      = cols_to_compare,
#         kind      = "scatter",
#         diag_kind = "hist",
#         height    = kwargs.get("height", 6),
#         aspect    = kwargs.get("height_multiplier", 1),
#         corner    = True,
#         dropna    = True,
#     )
    
#     if kwargs.get("write"):
#         output_image_dir = kwargs.get("output_image_dir", "for_paper_images")
#         plot_type        = kwargs.get("plot_type")
#         x                = kwargs.get("x")
#         runname          = kwargs.get("runname")
#         filetype         = kwargs.get("filetype", "png")
        
#         fig.savefig(f"{output_image_dir}/{plot_type}_{x}_{runname}.{filetype}")
        
#     if kwargs.get("show"):
#         fig.show()
        
    fig = px.scatter_matrix(
        df, 
        dimensions = cols_to_compare,
        color      = dict_o_kwargs.get("color"),
        labels     = labels, # remove underscore
        color_discrete_sequence = dict_o_kwargs.get("color_discrete_sequence"),
        color_discrete_map      = dict_o_kwargs.get("color_discrete_map"),
        opacity    = 0.4
    )
    
    fig.update_traces(
        diagonal_visible = False,
        #showupperhalf    = False,
    )

    return fig

#=========================================================================
#=========================================================================

def create_histogram(
    x, 
    df, 
    dict_o_kwargs, 
):
    
    df[x] +=  dict_o_kwargs.get("hist_offset", 0)
    fig = px.histogram(
        df,
        x                       = x,
        color                   = dict_o_kwargs.get("color"),
        color_discrete_sequence = dict_o_kwargs.get("color_discrete_sequence"),
        color_discrete_map      = dict_o_kwargs.get("color_discrete_map"),
        histnorm                = dict_o_kwargs.get("histnorm"),
        facet_col               = dict_o_kwargs.get("facet_col"),
        facet_row               = dict_o_kwargs.get("facet_row"),
        nbins                   = dict_o_kwargs.get("nbins", 0),
        marginal                = dict_o_kwargs.get("marginal"),
        #hover_data = {'Galaxy ID': (":c", full_df.index)},
    )
    
    if dict_o_kwargs.get("facet_col") or dict_o_kwargs.get("facet_row"):
        fig.for_each_annotation(lambda a: a.update(text = a.text.split("=")[-1]))

    # if multi:
    #     fig.update_layout(barmode = "overlay")
    #     fig.update_traces(
    #         opacity = 0.75,
    #         marker_line_width = 1,
    #         marker_line_color = "white"
    #     )
    
    return fig

#=========================================================================
#=========================================================================

def create_overlay_histogram(
    x,
    dfs,
    dict_o_kwargs
):
    
    histnorm    = dict_o_kwargs.get("histnorm", "probability")
    nbinsx      = dict_o_kwargs.get("nbins", 40)
    
    xmin, xmax  = None, None
    xaxis_range = dict_o_kwargs.get("xaxis_range", (xmin, xmax))
    
    if isinstance(xaxis_range, (tuple, list)) and all(xaxis_range):
        xmin, xmax  = xaxis_range
        
    xbins       = None
    
    # Because xmin could be 0
    if nbinsx and isinstance(xmin, (int, float)) and isinstance(xmax, (int, float)):
        #print(xmin, xmax, nbinsx)
        xbins = dict(
            start = xmin,
            end   = xmax,
            size  = (xmax - xmin)/nbinsx,
        )
    
    colors = deepcopy(px.colors.qualitative.Plotly)
    # Bulge -- Redder
    colors[0] = px.colors.qualitative.Plotly[1]
    # Disk -- Bluer
    colors[1] = px.colors.qualitative.Plotly[0]
    
    num_rows, num_cols = ceil(len(dfs)/3), len(dfs)
    # facet_col = dict_o_kwargs.get("facet_col")
    # if facet_col:
    #     num_cols  = len(df[facet_col].unique())
    #     reversed_cols.pop(facet_col)
        
    # facet_row = dict_o_kwargs.get("facet_row")
    # if facet_row:
    #     num_rows  = len(df[facet_row].unique())
    #     reversed_cols.pop(facet_row)
    
    # colors and to_plot is already reversed so they're fine as-is
    fig = make_subplots(
        num_rows, 
        num_cols,
        subplot_titles = list(dfs.keys()),
        shared_yaxes = True,
    )
    
    dfs = list(dfs.values())
    showlegend = True
    for row in range(num_rows):
        for col in range(num_cols):
            df = dfs[3*row + col]
            
            reversed_cols   = list(reversed(deepcopy(df.columns)))
            reversed_colors = reversed(colors[:len(df.columns)])
    
            for col_name, color in zip(reversed_cols, reversed_colors):
                name     = col_name.split("_")
                name[-1] = str(int(name[-1]) + 1) # since we index at 0
                name     = " ".join(name)
                
                                    
                fig.add_trace(
                    go.Histogram(
                        x                 = df[col_name] + dict_o_kwargs.get("hist_offset", 0),
                        histnorm          = histnorm,
                        name              = name,
                        marker_color      = color,
                        nbinsx            = nbinsx,
                        xbins             = xbins,
                        showlegend        = showlegend,
                        bingroup          = 1
                        #marginal          = dict_o_kwargs.get("marginal"),
                        #hover_data = {'Galaxy ID': (":c", full_df.index)},
                    ),
                    row = row + 1, col = col + 1
                )
            # Turning this off after one go-round
            showlegend = False
        
    # if facet_col or facet_row:
    #     fig.for_each_annotation(lambda a: a.update(text = a.text.split("=")[-1]))
        
    # This applies the yaxis title text to just 'one' of the facet cols, i.e. the first
    fig.update_layout(
        barmode           = "overlay",
        legend_traceorder = "reversed",
        yaxis_title_text  = histnorm,
    )
    
    fig.update_traces(
        opacity           = 0.75,
        marker_line_width = 1,
        marker_line_color = "white",
    )
    
    # This applies x to *all* facet cols
    fig.update_xaxes(
        title_text  = x,
        #minor_ticks = "outside",
    )
    
    # fig.update_yaxes(
    #     minor_ticks = "outside",
    # )
    
    return fig

#=========================================================================
#=========================================================================

def create_plot(
    x, 
    runname, 
    plot_type, 
    df_or_dfs, 
    output_image_dir = "for_paper_images", 
    **kwargs
):
        
    dict_o_kwargs = {
        "y"               : None,
        
        # Need these for histogram binning
        "xaxis_range"     : (None, None),
        "log_x"           : None,
        
        "color"                     : None,
        "color_order"               : None,
        "color_discrete_sequence"   : None,
        "color_discrete_map"        : None,
        "range_color"               : None,
        #"color_continuous_midpoint" : None,
        
        "hist_offset"     : 0,
        "histnorm"        : "",
        "marginal"        : None,
        "nbins"           : 0,
        
        "facet_col"       : None,
        "facet_row"       : None,
        
        "cutoff_val"      : 0.007,
        "add_vline"       : True,
        "add_hline"       : True,
    }
    
    # Updating with kwargs
    dict_o_kwargs = {key : kwargs.get(key, default) for key, default in dict_o_kwargs.items()}
    
    #plt.clf()
    #pio.templates.default = "plotly_white"
    
    if plot_type == "ecdf":
        fig = create_ecdf(x, df_or_dfs, dict_o_kwargs)
    elif plot_type == "scatter":
        fig = create_scatter(x, dict_o_kwargs.get("y"), df_or_dfs, dict_o_kwargs)
    elif plot_type == "scatter_matrix":
        fig = create_scatter_matrix(df_or_dfs, kwargs.get("cols_to_compare"), dict_o_kwargs)#, **kwargs)
        #return fig
    elif plot_type == "histogram":
        fig = create_histogram(x, df_or_dfs, dict_o_kwargs)
        #kwargs["yaxis_title"] = kwargs.get("yaxis_title", kwargs.get("histnorm", "probability"))
    elif plot_type == "overlay_histogram":
        fig = create_overlay_histogram(x, df_or_dfs, dict_o_kwargs)
    else:
        return
    
    if kwargs.get("title"):
        fig.update_layout(
            title_text = kwargs.get("title"), 
            title_x    = kwargs.get("title_x"), 
            title_y    = kwargs.get("title_y")
        )
    
    row, col = 1, None    
#     if plot_type != "overlay_histogram":
#         row, col = 1, 1

    # This applies x to *all* facet cols
    if "_" in x:
        fig.update_xaxes(
            title_text = " ".join(x.split("_")), row = row, col = col
        )
    
    # Only have one axis title instead of multiple stacked since this updates all cols as well
    y = dict_o_kwargs.get("y")
    if y and "_" in y:
        fig.update_yaxes(
            title_text = " ".join(y.split("_")), row = row, col = 1
        )
        
    if kwargs.get("xaxis_title"):
        fig.update_xaxes(title_text = kwargs.get("xaxis_title"), row = row, col = col)
    if kwargs.get("yaxis_title"):
        fig.update_yaxes(title_text = kwargs.get("yaxis_title"), row = row, col = col)
        
    if kwargs.get("log_x"):
        fig.update_xaxes(type = "log", row = row, col = col)
    if kwargs.get("log_y"):
        fig.update_yaxes(type = "log", row = row, col = col)
       
    xaxis_range = kwargs.get("xaxis_range")
    yaxis_range = kwargs.get("yaxis_range")
    if isinstance(xaxis_range, (tuple, list)):
        fig.update_xaxes(range = kwargs.get("xaxis_range"), row = row, col = col)
    if isinstance(yaxis_range, (tuple, list)):
        fig.update_yaxes(range = kwargs.get("yaxis_range"), row = row, col = col)
        
    fig.update_xaxes(minor_ticks = "outside", row = row, col = col)
    fig.update_yaxes(minor_ticks = "outside", row = row, col = col)
    
    height            = kwargs.get("height", 800)
    height_multiplier = kwargs.get("height_multiplier", 1.5) #1200
    width             = height*height_multiplier
    
    if kwargs.get("show"):
        # Invert boolean so interactive = True means NOT static plot and vice-versa
        fig.show(
            config = {
                'staticPlot' : not kwargs.get("interactive", False),
                'toImageButtonOptions': {
                    'height': height,
                    'width' : width
                }
            }
        )
    
    if kwargs.get("write"):
        filetype = kwargs.get("filetype", "png")
        fig.write_image(
            f"{output_image_dir}/{plot_type}_{x}_{runname}.{filetype}", 
            height = height, 
            width = width
        )
        
    fig.data   = []
    fig.layout = {}
    
    return fig
