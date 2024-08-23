from pandas import read_csv
from os.path import exists
from numpy import float32

def load_petromags(
    gzoo_file,
    color_band
):
    
    if not exists(gzoo_file):
        return full_df
    
    petromag_col = f"petroMag_{color_band}"
    gname_col = "GZ_dr8objid"
    
    gzoo_data = read_csv(
        gzoo_file, 
        sep = "\t", 
        index_col = gname_col,
        usecols  = [gname_col, petromag_col],
        dtype = {
            gname_col    : str,
            petromag_col : float32
        }
    ).fillna("None")
    
    gzoo_data = gzoo_data[~gzoo_data.index.duplicated(keep='first')]
    
    return gzoo_data
