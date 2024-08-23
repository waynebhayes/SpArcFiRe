from astropy.io import fits

def grab_header_parameters_for_flux(input_file):
    try:
        with fits.open(input_file) as hdul:
            exptime = hdul[1].header.get("EXPTIME", 1)
            mag_zpt = hdul[2].header.get("MAGZPT", 24.8) # Default for SDSS r band
            
    except FileNotFoundError as fe:
        print(f"Could not fits.open {input_file}.")
        print("Using exptime and mag zpt defaults from SDSS DR7 r-band.")
        exptime, mag_zpt = 1, 24.8
        
    return exptime, mag_zpt