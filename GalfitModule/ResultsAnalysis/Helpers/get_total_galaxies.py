from Functions.helper_functions import find_files

def get_total_galaxies(in_dir = "sparcfire-in", out_dir = "sparcfire-out"):   
    all_gnames_in  = find_files(in_dir, "123*", "f")
    all_gnames_out = find_files(out_dir, "123*", "d")
    total_galaxies = min(len(all_gnames_in), len(all_gnames_out))
    if not total_galaxies:
        total_galaxies  = max(len(all_gnames_in), len(all_gnames_out))
        
    return total_galaxies