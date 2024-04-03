import shutil
import sys
import os
from os.path import join as pj
# from os.path import exists
from astropy.io import fits
from collections import ChainMap

#from joblib import cpu_count #Parallel, delayed
import asyncio
import psutil
from itertools import permutations
#from Collections import namedtuple as nt

_HOME_DIR = os.path.expanduser("~")    

try:
    _SPARCFIRE_DIR = os.environ["SPARCFIRE_HOME"]
    _MODULE_DIR = pj(_SPARCFIRE_DIR, "GalfitModule")
except KeyError:
    # print("SPARCFIRE_HOME is not set. Please run 'setup.bash' inside SpArcFiRe directory if not done so already.")
    # print("Running on the assumption that GalfitModule is in your home directory... (if not this will fail and quit!)") 
    _MODULE_DIR = pj(_HOME_DIR, "GalfitModule")
    
sys.path.append(_MODULE_DIR)
    
#from Classes.Components import *
#from Classes.Containers import *
# FitsHandlers imports all of the above
from sparc_to_galfit_feedme_gen import *
from Classes.FitsHandlers import *
from Functions.helper_functions import *
#from Utilities.parallel_residual_calc import fill_objects #, parallel_wrapper

import star_removal.no_log_remove_stars_with_sextractor as remove_stars_with_sextractor

# ns2_component_switches = nt("component_switches", ["bulge", "disk_for_arms", "arms", "fourier", "sky"])

# ns3_component_switches = nt("component_switches", ["bulge", "disk", "disk_for_arms", "arms", "fourier", "sky"])

def rerun_galfit(galfit_output, base_galfit_cmd, *args):
    
    if galfit_output.success:   
        # don't use this anymore
        #galfit_num = f"{int(galfit_output.galfit_num) + 1:0>2}"
        
        print("Re-running GALFIT...")
        # String comparison until I implement subtraction ;)
        # Use this to generically determine which components were optimized
        updated_components = (comp for comp, icomp in zip(galfit_output.to_tuple(), args) if str(comp) != str(icomp))

        # run_galfit_cmd = f"{base_galfit_cmd} galfit.{galfit_num}"
        galfit_output.to_file(*updated_components)
        
        run_galfit_cmd = f"{base_galfit_cmd} {galfit_output.path_to_feedme}"

        #galfit_num = f"{int(galfit_output.galfit_num) + 1:0>2}"
        galfit_output = OutputContainer(sp(run_galfit_cmd, capture_output = True), **galfit_output.to_dict())
       
    else:
        print("Galfit failed!")

    return galfit_output

def touch_failed_fits(gname, tmp_fits_dir):
    
    prefix = "failed_"
    #filepath = os.path.dirname(tmp_fits_dir)
    fail_fileout = pj(tmp_fits_dir, f"{prefix}{gname}_galfit_out.fits")
    #sp(f"touch {fail_fileout}", capture_output = False)
    with open(fail_fileout, mode='a'): pass
    
def check_for_starmasks(gnames, masks_dir):
    masks_in_dir = [
        os.path.basename(i).split("_star-rm.fits")[0] 
        for i in find_files(masks_dir, "*star-rm.fits", "f")
    ]
    return list(set(gnames).difference(set(masks_in_dir)))

def check_success(in_fits, previous_state = False):
    if not isinstance(in_fits, OutputContainer):
        print("Not an OutputFits type!")
    
    else:
        if in_fits.success:
            return True
    
    return (previous_state or False)

def generate_model_for_sparcfire(
    gfits, 
    base_galfit_cmd, 
    in_dir, 
    tmp_fits_dir = "",
    output_dir   = "",
    rad          = (None, None)
):
    gname  = os.path.basename(gfits).replace("_galfit_out.fits", "")
    suffix = "for_sparcfire"
    if not tmp_fits_dir:
        tmp_fits_dir = os.path.split(gfits)[0]
        
    fits_file = OutputFits(gfits, load_default = False)
    feedme = fits_file.feedme
    
    feedme.header.input_image.value   = pj(in_dir, f"{gname}.fits")
    feedme.header.output_image.value  = pj(tmp_fits_dir, f"{gname}_{suffix}.fits")
    if output_dir:
        feedme.header.output_image.value  = pj(output_dir, f"{gname}_{suffix}.fits")

    center = feedme.sersic_0.position
    # In an attempt to avoid more I/O...
    # Alas I don't think I can use this
    rad_x  = rad[0]
    rad_y  = rad[1]
    
    if not all(rad):
        # No crop to allow SpArcFiRe to use its autocrop routine 
        # without issue
        with fits.open(pj(in_dir, f"{gname}.fits")) as in_fits:
            input_file_shape = np.shape(in_fits[0].data)
            rad_x, rad_y = input_file_shape[0] // 2, input_file_shape[1] // 2

    feedme.header.region_to_fit.value = (
        center.x - rad_x, 
        center.x + rad_x, 
        center.y - rad_y, 
        center.y + rad_y
    )
    feedme.header.optimize.value      = 1
    
    if "power_0" in feedme.components: #or "arms" in feedme.components:
        power_comp_number = feedme.power_0.component_number
        # Strongly detuning the arms to give SpArcFiRe the best chance of success
        rot_comp = [
            name for name, comp in feedme.components.items()
            if  comp.component_number == power_comp_number
            and comp.component_type   == "sersic"
        ][0]
        
        feedme.components[rot_comp].axis_ratio.value = 0.15
        
        # Set all other sersic components to 0 to just keep the disk being rotated
        # and the header and sky
        for c_name, comp in feedme.components.items():
            if "sky" in c_name or "header" in c_name:
                continue
                
            if comp.component_number != power_comp_number:
                comp.skip.value = 1
        
    filename = pj(tmp_fits_dir, f"{gname}_{suffix}.in")
    feedme.to_file(filename = filename)
    _ = sp(f"{base_galfit_cmd} {filename}", capture_output = True)
    
    return feedme.header.output_image.value

async def async_sp(*args, **kwargs):
    
    # BIG thanks to these blogs
    # https://blog.dalibo.com/2022/09/12/monitoring-python-subprocesses.html
    # https://blog.meadsteve.dev/programming/2020/02/23/monitoring-async-python/
    
    timeout = kwargs.pop("timeout", 3*60)
    
    # Some defaults in case I've implemented things incorrectly...
    stdout = "...now exiting to system..."
    stderr = ""
    
    #proc = await asyncio.create_subprocess_shell(
    #args = [a.split() for a in args][0]
    proc = await asyncio.create_subprocess_exec(
        *[a.split() for a in args][0],
        stdout = asyncio.subprocess.PIPE,
        stderr = asyncio.subprocess.PIPE,
        **kwargs
    )
    
    async def check_time():
        loop = asyncio.get_running_loop()
        while loop.is_running():

            # This can happen if it's running too fast ;) 
            # i.e. for some Bulge fits
            try:
                #print(f"Time {psutil.Process(proc.pid).cpu_times().user}")
                exec_time = psutil.Process(proc.pid).cpu_times().user
            except psutil.NoSuchProcess:
                break
            
            if exec_time >= timeout:
                print(f"Subprocess call timed out by cpu time.")
                #stdout = "...now exiting to system..."
                stderr = f"GALFIT timed out."
                #raise TimeoutError(f"GALFIT timed out.")
                proc.kill()
                #await proc.wait()
                return stdout, stderr
            
            # Check every 30 seconds
            await asyncio.sleep(30)
            
    #try:
        #stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout = timeout)
        # p = psutil.Process(proc.pid)
        # user_time = p.cpu_times().user
        
    ct = asyncio.create_task(check_time())
    #(stdout, stderr), _ = await asyncio.gather(proc.communicate(), ct, return_exceptions = True)
    (stdout, stderr), timeout_tup = await asyncio.gather(proc.communicate(), ct, return_exceptions = True)

    # This can only be used for python>3.11 and theoretically works better than the current method
    #async with asyncio.timeout(timeout):
    #    stdout, stderr = await proc.communicate()

    if timeout_tup:
        stdout = timeout_tup[0]
        stderr = timeout_tup[1]

    if isinstance(stdout, bytes):
        stdout = stdout.decode('utf8')

    if isinstance(stderr, bytes):
        stderr = stderr.decode('utf8')
        
#     except TimeoutError:
#         #if proc.returncode is None:
#         #    parent = psutil.Process(proc.pid)
#         #    for child in parent.children(recursive=True): 
#         #        child.terminate()
#         #    parent.terminate()

#         #await asyncio.sleep(1)
#         proc.kill()
#         await proc.wait()
#         #process._transport.close()
#         print(f"Subprocess call timed out.") 
#         #print(f"Timeout is set to {timeout/60} minutes. If this is not enough, please consider increasing this value.")
#         #if kwargs.get("verbose"):
        
#         stdout = "...now exiting to system..."
#         stderr = ""
#         #reset_all()

    # except Exception as e:
    #     print(f"Something went wrong in grabbing stdout from GALFIT! {e}")
    #     stdout = "...now exiting to system..."
    #     stderr = ""
    
    # For debuggin
    #print(*stdout.decode('utf8').split("\n")[:75], sep = "\n")
    
    return subprocess.CompletedProcess("", 0, stdout = stdout, stderr = stderr)

async def run_galfit(
    component_switches,
    galfit_output, 
    tmp_feedme_in, 
    base_galfit_cmd, 
    **kwargs
):
    
    timeout      = kwargs.get("timeout", 3*60)
    # To be abundantly safe, this is turned off
    load_default = kwargs.get("load_default", False)
    
    header = galfit_output.header
    
    # Stick to the convention of assuming fourier and sky
    # switch[0] is the switch to use the component
    filepath_in = tmp_feedme_in
    if kwargs.get("add_suffix", True):
        suffix = "+".join([
            cname for cname, switch in component_switches.items() 
            if switch[0] and cname not in ("header", "fourier", "sky")
        ])
        filepath_in = f"{tmp_feedme_in.replace('.in', '')}_{suffix}.in"
    
    # Equivalent to initial_components/galfit_output.components[cname]
    # switch[1] is the switch to use initial (0) or galfit_output (1)
    comps_to_file = {
        cname : galfit_output.components[cname] for cname, switch in component_switches.items() if switch[0]
    }
    
    # to_file no arguments uses the feedme path attribute for location
    # and uses all the components in the container object (in whatever order they're in)
    header.to_file(filepath_in, *comps_to_file.values())
    run_galfit_cmd = f"{base_galfit_cmd} {filepath_in}"

    # We can use galfit_output here regardless since we have deepcopied it
    # from initial components
    if kwargs.get("use_async", True):
        galfit_output = OutputContainer(
            await async_sp(run_galfit_cmd, timeout = timeout),
            path_to_feedme = filepath_in,
            load_default   = load_default,
            store_text     = True,
            header         = galfit_output.header,
            **comps_to_file #galfit_output.components
        )

    else:
        galfit_output = OutputContainer(
            sp(run_galfit_cmd), #, timeout = timeout), 
            path_to_feedme = filepath_in,
            load_default   = load_default,
            store_text     = True,
            header         = galfit_output.header,
            **comps_to_file #galfit_output.components
        )
        
    # Turn on the switches for using galfit_output if run this round
    component_switches.update({
        cname : [1, 1] for cname, switch in component_switches.items() if switch[0]
    })
    
    if kwargs.get("verbose", False):
        print(galfit_output)#.components.values())

    return galfit_output, component_switches
    
async def parameter_search_fit(
    bulge_magnitude, # loop variables
    disk_magnitude,  # loop variables
    #bulge_sersic,
    gname,
    initial_feedme, 
    base_galfit_cmd,
    disk_axis_ratio,
    **kwargs
):
    # With this timeout set, galaxies will take no longer than two hours per
    # for default settings searching 20 parameters synchronously
    # i.e. 2 steps, 20 parameter variations, 1 minutes max per = 40 min.
    
    # The value timeout is set to should, of course, depend on the average 
    # image resolution and expected runtime
    timeout   = 60 # Seconds
    timeout  *= 3  # Minutes
    #psutil.Process(proc.pid).cpu_times().user >= timeout:
    
    out_dir        = kwargs.get("out_dir")
    num_steps      = int(kwargs.get("num_steps", 2))
    num_components = int(kwargs.get("num_components", 3))
    
    #tmp_fits_dir  = kwargs.get("tmp_fits_dir")
    
    # This doesn't change
    # initial_feedme is a FeedmeContainer object
    #header       = initial_feedme.header
    #feedme_path  = initial_feedme.path_to_feedme
    output_image = initial_feedme.header.output_image.value
    
    # This is to clean up the code and make it easier to understand
    # Make a deepcopy just in case
    initial_components = deepcopy(initial_feedme)

    use_spiral = True
    if "arms" not in initial_components.components and "power" not in " ".join(initial_components.components.keys()):
        use_spiral = False
    
    header = initial_components.header
    
    initial_components.bulge.magnitude.value += bulge_magnitude
    #initial_components.bulge.magnitude.value = bulge_magnitude
    #initial_components.disk.magnitude.value  = disk_magnitude
    
    #if use_spiral:
    #    initial_components.disk_for_arms.magnitude.value  = disk_magnitude

    # ---------------------------------------------------------------------
    # Note to self, OutputContainer is *just* for handling output
    # It does not retain the input state fed into Galfit
    # And the input dict is *before* output but will be updated *by* output
    # ---------------------------------------------------------------------
    
    bulge_str = f"m{10*bulge_magnitude:.0f}"
    #disk_str  = f"m{disk_magnitude}"
    #disk_str  = "m" + str(bulge_sersic).replace(".", "")
                
    #print(f"Bulge magnitude {bulge_magnitude}, Disk magnitude {disk_magnitude}")
    
    basepath, basename = os.path.split(output_image)
    #new_basename       = f"{bulge_str}{disk_str}_{basename}"
    new_basename       = f"{bulge_str}_{basename}"
    new_output_image   = pj(basepath, new_basename)
    initial_components.header.output_image.value = new_output_image
    
    feedme_name     = os.path.basename(initial_feedme.path_to_feedme)
    #new_feedme_name = f"{bulge_str}{disk_str}_{feedme_name}"
    new_feedme_name = f"{bulge_str}_{feedme_name}"
    tmp_feedme_in   = pj(basepath, new_feedme_name)

    success = False
    
    # This needs to be done here so that we can do a proper 1/2/3 step fit
    galfit_output = deepcopy(initial_components)
    
    # first number is switch for using bulge, disk_for_arms, arms, fourier, sky
    # second is switch for taking from original feedme or previous output... TODO: remove?
    # Starting choices
    if num_steps == 1:
        component_switches = {cname : [1,0] for cname in initial_components.components.keys()}
    else:
        component_switches = {cname : [0,0] for cname in initial_components.components.keys()}
        
        component_switches["bulge"] = [1,0]
        component_switches["sky"]   = [1,0]
            
    # bulge, disk, disk_for_arms, arms, fourier, sky
    if num_components == 3 and num_steps == 2:
        # Most of these are turned on, might as well re-dict-comp
        component_switches = {cname : [1,0] for cname in initial_components.components.keys()}
        component_switches["disk"] = [0,0]
            
    component_switches.pop("header")
    
    # Any galfit runs that feed into an output container *must* have capture_output = True or 
    # they won't update the new parameters
    # TODO: Add timeout to all galfit sp calls (5 min) (???)
    
    # ns2, nc2, spiral ==    bulge, bulge+disk+arms
    # ns3, nc2, spiral ==    bulge, bulge+disk, bulge+disk+arms
    #
    # ns2, nc2, no spiral == bulge, bulge+disk (REPEAT)
    # ns3, nc2, no spiral == (REPEAT)
    #
    # ns2, nc3, spiral ==    bulge+arm_disk+arms, bulge+disk+arm_disk+arms
    # ns3, nc3, spiral ==    bulge, bulge+arm_disk+disk, bulge+disk+arm_disk+arms
    #
    # ns2, nc3, no spiral == (REPEAT)
    # ns3, nc3, no spiral == (REPEAT)
    
    if not use_spiral:
        if num_steps >= 2:
            # STEP 1 / 2 or 3
            galfit_output, component_switches = await run_galfit(
                component_switches,
                galfit_output, 
                tmp_feedme_in, 
                base_galfit_cmd, 
                **kwargs
            )
            success = check_success(galfit_output, success)
            
            galfit_output.add_components(disk = initial_components.disk)
    
    # bulge -> bulge + disk + arms
    elif num_steps >= 2:
        
        # Number of steps 2: bulge -> bulge+disk_for_arms+arms
        if num_components == 2:
            
            # Only turn load default on if order doesn't matter
            # In the case of num_components = 3, it *does* matter
            # Since we assume continugous component nubmer when 
            # extracting from output
            
            # STEP 1 / 2 or 3
            galfit_output, component_switches = await run_galfit(
                component_switches,
                galfit_output, 
                tmp_feedme_in, 
                base_galfit_cmd,
                #load_default = True,
                **kwargs
            )
            success = check_success(galfit_output, success)
            
            galfit_output.add_components(disk_for_arms = initial_components.disk_for_arms)
            # Number of steps 3: bulge -> bulge+disk_for_arms -> bulge+disk_for_arms+arms
            
            # STEP 2 / 3
            if num_steps == 3:
                # Use initial feedme since it hasn't been used yet
                component_switches["disk_for_arms"] = [1,0]
                
                galfit_output, component_switches = await run_galfit(
                    component_switches,
                    galfit_output, 
                    tmp_feedme_in, 
                    base_galfit_cmd,
                    load_default = True,
                    **kwargs
                )
                
                success = check_success(galfit_output, success)
                
            # Prep for spiraling
            galfit_output.disk_for_arms.axis_ratio.value = disk_axis_ratio
            galfit_output.add_components(arms = initial_components.arms, fourier = initial_components.fourier)
                
        # Number of steps 2: bulge+disk_for_arms+arms -> bulge+disk+disk_for_arms+arms
        elif num_components == 3:
            
            galfit_output.disk_for_arms.axis_ratio.value = disk_axis_ratio
            
            if num_steps == 2:
                
                # STEP 1 / 2
                galfit_output, component_switches = await run_galfit(
                    component_switches,
                    galfit_output, 
                    tmp_feedme_in, 
                    base_galfit_cmd, 
                    **kwargs
                )
                success = check_success(galfit_output, success)
                
            # Number of steps 3: bulge -> bulge+disk_for_arms -> bulge+disk_for_arms+arms
            elif num_steps == 3:
                
                # STEP 1 / 3
                galfit_output, component_switches = await run_galfit(
                    component_switches,
                    galfit_output, 
                    tmp_feedme_in, 
                    base_galfit_cmd, 
                    **kwargs
                )
                
                success = check_success(galfit_output, success)
                
                # STEP 2 / 3
                # Use initial feedme since it hasn't been used yet
                component_switches["disk_for_arms"] = [1,0]
                component_switches["arms"]          = [1,0]
                component_switches["fourier"]       = [1,0]
                
                galfit_output, component_switches = await run_galfit(
                    component_switches,
                    galfit_output, 
                    tmp_feedme_in, 
                    base_galfit_cmd, 
                    **kwargs
                )
                
                success = check_success(galfit_output, success)

            galfit_output.add_components(disk = initial_components.disk)
            galfit_output.disk.magnitude.value = galfit_output.disk_for_arms.magnitude.value 
        
        # Only fix sky if first step is successful
        # if galfit_output.success:
        #     # Fix sky parameters per Galfit 'tips' recommendation
        #     for key in galfit_output.sky.param_fix:
        #         galfit_output.sky.param_fix[key] = 0
        
            # After determining sky and initial component(s) with extended background,
            # shrink fitting region closer to the galaxy itself
            # xcenter, ycenter = galfit_output.bulge.position
            # old_region = galfit_output.header.region_to_fit
            # crop_mult = 2
            # crop_rad = 1.5*(xcenter - old_region[0])/crop_mult
            # galfit_output.header.region_to_fit = (
            #                               max(round(xcenter - crop_rad), 0), # Just in case
            #                               round(xcenter + crop_rad), 
            #                               max(round(ycenter - crop_rad), 0), 
            #                               round(ycenter + crop_rad)
            #                              )
            # galfit_output.header.update_param_values()

            #galfit_output.arms.param_fix["outer_rad"] = 1
    
    # FINAL OUTPUT
    # Turn on all switches
    # Using switch[1] is redundant but again, just in case
    component_switches = {cname : [1, switch[1]] for cname, switch in component_switches.items()}
    
    if use_spiral and num_steps == 1:
        galfit_output.disk_for_arms.axis_ratio.value = disk_axis_ratio
    
    galfit_output, component_switches = await run_galfit(
        component_switches,
        galfit_output, 
        tmp_feedme_in, 
        base_galfit_cmd,
        add_suffix = False,
        **kwargs
    )
                
    success = check_success(galfit_output, success)
    
    if kwargs.get("verbose") and not use_spiral:
        print(str(galfit_output))

    # ********************************************************************
    # Release holds on inner and outer spiral radius
    # does not seem to make a significant difference
    # ********************************************************************
    if use_spiral:
        galfit_output.arms.inner_rad.fix = 1
        galfit_output.arms.outer_rad.fix = 1

        final_galfit_output, component_switches = await run_galfit(
            component_switches,
            galfit_output, 
            tmp_feedme_in, 
            base_galfit_cmd,
            add_suffix = False,
            **kwargs
        )
        
        success = check_success(final_galfit_output, success)
        
        if kwargs.get("verbose"):
            print(str(final_galfit_output))
        
    if success:
        
        fits_file = OutputFits(new_output_image)
        use_bulge_mask = False #True
        
        try:
            mask_fits_file = FitsFile(header.pixel_mask.value)
        except Exception as e:
            print(f"There was an issue opening star mask for {new_basename}. Proceeding without mask...")
            # Logic implemented to handle None
            mask_fits_file = None #np.zeros((500,500))

        # If skip is enabled then arms are not fit so bulge masking doesn't make sense
        c_types = [comp.component_type for comp in fits_file.feedme.components.values()]
        
        if use_bulge_mask and out_dir and ("power" in c_types):
            _ = fits_file.generate_bulge_mask(pj(out_dir, gname, f"{gname}.csv"))
        else:
            use_bulge_mask = False

        masked_residual_normalized = fits_file.generate_masked_residual(
            mask_fits_file, 
            use_bulge_mask     = use_bulge_mask,
            update_fits_header = False
        )
        
        if masked_residual_normalized is None:
            print(f"Could not calculate nmr for {new_basename}.")
            # returning some absurd number
            return {new_output_image : {"pvalue" : 0, "nmr" : 100000}}

        # output image path : {p value : KS test p value, nmr : NMR})
        return {new_output_image : {"pvalue" : fits_file.kstest.pvalue, "nmr" : fits_file.nmr}}
    
    else:
        print(f"Could not calculate nmr for {new_basename}.")
        
    return {new_output_image : {"pvalue" : 0, "nmr" : 100000}}
    
async def wrapper(
    b_d_magnitudes,
    #b_magnitudes,
    gname,
    initial_feedme, 
    base_galfit_cmd,
    disk_axis_ratio,
    **kwargs
):
    
    # Even if we are not using async, I believe we must still "await" or an error will be thrown
    #loop = asyncio.get_event_loop()
    
    fitted_galaxies = await asyncio.gather(*(
        parameter_search_fit(
            bulge_magnitude,
            disk_magnitude,
            #bulge_sersic,
            gname,
            initial_feedme,
            base_galfit_cmd,
            disk_axis_ratio,
            **kwargs
        ) for (bulge_magnitude, disk_magnitude) in b_d_magnitudes
                          ), return_exceptions = True)
    
    # Convert list of dicts to single dict
    # gpath : nmr value
    try:
        fitted_galaxies = dict(ChainMap(*fitted_galaxies[::-1]))
    except TypeError as te: #(IndexError, AttributeError):
        if "IndexError" in str(te):
            fitted_galaxies = {k : v for k, v in fitted_galaxies.items() 
                               if v not in (IndexError, AttributeError)
                              }
        if fitted_galaxies:
            fitted_galaxies = dict(ChainMap(*fitted_galaxies[::-1]))
        else:
            print(f"Something went wrong, {gname} errored out. Counting as a failure.")
        
    return fitted_galaxies
    
def main(**kwargs):
    
    # Main Directories
    cwd       = kwargs.get("cwd"    , os.getcwd())
    in_dir    = kwargs.get("in_dir" , pj(cwd, "sparcfire-in"))
    tmp_dir   = kwargs.get("tmp_dir", pj(cwd, "sparcfire-tmp"))
    out_dir   = kwargs.get("out_dir", pj(cwd, "sparcfire-out"))
    
    # Should some directory structure change in the future
    tmp_fits_dir  = kwargs.get("tmp_fits_dir" , pj(tmp_dir, "galfits"))
    tmp_png_dir   = kwargs.get("tmp_png_dir"  , pj(tmp_dir, "galfit_png"))
    tmp_masks_dir = kwargs.get("tmp_masks_dir", pj(tmp_dir, "galfit_masks"))
    tmp_psf_dir   = kwargs.get("tmp_psf_dir"  , pj(tmp_dir, "psf_files"))
    
    out_png_dir   = kwargs.get("out_png_dir"  , pj(out_dir, "galfit_png"))
    
    # Of course the important things
    num_steps      = int(kwargs.get("num_steps", 2))
    num_components = int(kwargs.get("num_components", 3))
    
    parallel  = int(kwargs.get("parallel", 1))
    use_async = kwargs.get("use_async", True)
    
    # For simultaneous fitting
    simultaneous_fitting = kwargs.get("simultaneous_fitting", False)
    sim_fitting_dir  = kwargs.get("sim_fitting_dir", pj(tmp_dir, "sim_fitting"))
    
    sf_in_dir    = kwargs.get("sim_fitting_in_dir", pj(sim_fitting_dir, "sparcfire-in"))
    sf_tmp_dir   = kwargs.get("sim_fitting_tmp_dir", pj(sim_fitting_dir, "sparcfire-tmp"))
    sf_out_dir   = kwargs.get("sim_fitting_out_dir", pj(sim_fitting_dir, "sparcfire-out"))
    sf_masks_dir = kwargs.get("sim_fitting_masks_dir", pj(sf_tmp_dir, "galfit_masks"))
    
    # For verbosity, default to capturing output
    # Keep both for clarity
    verbose        = kwargs.get("verbose", False)
    capture_output = kwargs.get("capture_output", True)
    
    # For generating starmasks per galaxy
    generate_starmasks = kwargs.get("generate_starmasks", True)
    
    # To indicate running from a local tmp dir
    run_from_tmp       = kwargs.get("run_from_tmp", False)
    
    # To aggressively free up space
    aggressive_clean   = kwargs.get("aggressive_clean", False)
    
    # Feeding in as comma separated galaxies
    # If single galaxy, this returns a list containing just the one galaxy
    galaxy_names = kwargs.get("galaxy_names", [])
    # petromags         = kwargs.get("petromags", [])
    # bulge_axis_ratios = kwargs.get("bulge_axis_ratios", [])
    if isinstance(galaxy_names, str):
        galaxy_names = galaxy_names.split(",")
        
    if not galaxy_names:
        print("No galaxy names fed into go_go_galfit.")
        return None
        # petromags    = petromags.split(",")
        # bulge_axis_ratios    = bulge_axis_ratios.split(",")
    #     galaxy_names = [galaxy_names]
    #assert isinstance(galaxy_names, list), "input_filenames must be a list, even if it's a single galaxy."
    
    # This PITA brought to you by running things in the notebook environment
    run_galfit, run_fitspng, run_python = check_programs()
    run_galfit  = kwargs.get("run_galfit" , run_galfit)
    run_fitspng = kwargs.get("run_fitspng", run_fitspng)
    run_python  = kwargs.get("run_python" , run_python)
    
    # We could chunk this up for parallel but since Wayne's script 
    # automatically distributes processes, we don't have to worry about it
    # Deprecated
    # gname = ""
    # parallel = False
    # if len(galaxy_names) == 1:
    #     #gname = galaxy_names[0]
    #     parallel = True
    if run_from_tmp:
        for k,v in kwargs.items():
            if "_dir" in k and k != "out_dir":
                #_ = sp(f"rm -rf {v}")
            #os.makedirs(v, exist_ok = True)
                _ = sp(f"mkdir -p {v}")
    
    if generate_starmasks:
        print("Generating Starmasks (if need be)")
        
        galaxies_to_mask = check_for_starmasks(galaxy_names, tmp_masks_dir)
        
        if galaxies_to_mask:
            # I think we have to change path for source extractor
            star_removal_path = pj(_MODULE_DIR, "star_removal")
            os.chdir(star_removal_path)
            remove_stars_with_sextractor.main(in_dir, tmp_masks_dir, galaxies_to_mask)
            os.chdir(cwd)
        
    print("Running feedme generator...")
    feedme_info = write_to_feedmes(
        top_dir        = cwd,
        galaxy_names   = galaxy_names,
        in_dir         = in_dir,
        tmp_dir        = tmp_dir,
        out_dir        = out_dir,
        num_components = num_components
        # petromags = petromags,
        # bulge_axis_ratios = bulge_axis_ratios
    )
    
    if simultaneous_fitting:
        
        galaxies_to_mask = check_for_starmasks(galaxy_names, sf_masks_dir)
        
        if galaxies_to_mask:
            star_removal_path = pj(_MODULE_DIR, "star_removal")
            os.chdir(star_removal_path)
            remove_stars_with_sextractor.main(sf_in_dir, sf_masks_dir, galaxies_to_mask)
            os.chdir(cwd)
        
        sf_feedme_info = write_to_feedmes(top_dir = cwd,
                                          galaxy_names = galaxy_names,
                                          in_dir  = sf_in_dir,
                                          tmp_dir = sf_tmp_dir,
                                          out_dir = sf_out_dir,
                                          num_components = num_components
                                          # petromags = petromags,
                                          # bulge_axis_ratios = bulge_axis_ratios
                                          )
                                   
    max_it = 200
    base_galfit_cmd = f"{run_galfit} -imax {max_it}"
    # TODO: Move output sigma file as well, figure out race condition
    # base_galfit_cmd = f"{run_galfit} -imax {max_it} -outsig"
    
    failed          = []
    fitted_galaxies = []
    
    # Non-inclusive of end
    # I keep b_d_magnitudes even though I fix the disk for ease of modification
    # in the future
    b_d_magnitudes = [(b, 15) for b in np.linspace(-2, 2, num = 5)]
    # b_d_magnitudes = [(b, 15) for b in range(12, 17)]
    #b_magnitudes = range(12, 17)
    #b_d_magnitudes = [(b, d) for b in range(12, 17) for d in range(13, 17)]
    #[(b, d*0.1) for b in range(12, 17) for d in range(5, 41, 5)]
    
    # Working on non-parallel for now
    print()
    print(f"Galfitting with {num_components} components and {num_steps} component steps... (stdout is captured):")
    for gname in galaxy_names:
        # For debugging
        #feedme_info[gname].header.optimize.value = 1
        #if gname != "1237668298219127022":
        #    continue
        
        if gname not in feedme_info:
            print(f"Feedme not generated for {gname}")
            touch_failed_fits(gname, tmp_fits_dir)
            failed.append(gname)
            continue
            
        print(gname)
        
        disk_axis_ratio = 0.7
        with fits.open(pj(in_dir, f"{gname}.fits")) as gf:
            try:
                #SURVEY = SDSS-r  DR7
                #color = gf[0].header["SURVEY"].split()[0][-1]
                spirality_str = f"spirality"
                p_spirality = float(gf[0].header.get("spirality", 0))
                # p_spirality, disk_axis_ratio
                if p_spirality >= 0.9: 
                    disk_axis_ratio = 0.5
                elif p_spirality >= 0.75:
                    disk_axis_ratio = 0.6
                    
            except KeyError:
                print(f"There is no spirality keyword in the header for {gname}.")
        
        # ======================================== BEGIN GALFIT PARAMETER SEARCH LOOP ========================================    
        
        # FOR DEBUGGING 
        # kwargs["use_async"] = False
        # use_async = False
        # b_d_magnitudes = [(16,16)]
        #if parallel in (0, 1):
        #    use_async = True
        
        if use_async:
            # Limiting our # of asynchronous processes
            chunk = len(b_d_magnitudes)
            #chunk = len(b_magnitudes)
            if parallel == 1: #in (1, 2):
                chunk = 5

            #if len(b_d_magnitudes) > chunk:
            fitted_galaxies = {}
            for i in range(0, len(b_d_magnitudes), chunk):
                if (i + chunk) < len(b_d_magnitudes):
                    b_d_chunk = b_d_magnitudes[i:i + chunk]
                else:
                    b_d_chunk = b_d_magnitudes[i:]

                fitted_galaxies.update(asyncio.run(wrapper(
                    b_d_chunk,
                    gname,
                    feedme_info[gname], 
                    base_galfit_cmd,
                    disk_axis_ratio,
                    **kwargs
                )))
                
        else:
            fitted_galaxies = [asyncio.run(parameter_search_fit(
                bulge_magnitude,
                disk_magnitude,
                #bulge_sersic,
                gname,
                feedme_info[gname],
                base_galfit_cmd,
                disk_axis_ratio,
                **kwargs
            )) for (bulge_magnitude, disk_magnitude) in b_d_magnitudes]

            # Convert list of dicts to single dict
            # gpath : nmr value
            fitted_galaxies = dict(ChainMap(*fitted_galaxies[::-1]))
 
        # ======================================== END GALFIT MAGNITUDE LOOP ========================================
        # if not any(fitted_galaxies):
        #     print("No successful deconvolutions.")
        #     continue
            
        tmp_png_path  = pj(tmp_png_dir, gname)
        tmp_fits_path_gname = pj(tmp_fits_dir, f"{gname}_galfit_out.fits")
        
        # nmr_x_1-p
        fitted_galaxies_nmr_x_1_p = {
            gpath : (1 - inner_dict["pvalue"])*inner_dict["nmr"] 
            for gpath, inner_dict in fitted_galaxies.items()
        }
        #print(fitted_galaxies_nmr_x_1_p)
        
        try:
            best_fit_gpath = min(fitted_galaxies, key = fitted_galaxies_nmr_x_1_p.get)
            
        except ValueError as ve:
            #print(ve)
            print(f"Something went wrong with galaxy {gname}, can't find best fit from parameter search.")
        else:
            # If filling each of the FITS files is desired... it may not be if too much I/O is happening
            
            best_fit_p   = fitted_galaxies[best_fit_gpath]["pvalue"]
            best_fit_nmr = fitted_galaxies[best_fit_gpath]["nmr"]
            
            try:
                with fits.open(best_fit_gpath, mode = "update", output_verify = "ignore") as hdul:
                    hdul[2].header["NMR"]      = (round(best_fit_nmr, 8), "Norm of the masked residual")

                    # pvalue is sometimes none but round can't handle it
                    if isinstance(best_fit_p, float):
                        hdul[2].header["KS_P"] = (round(best_fit_p, 8), "p value of kstest vs noise")
                    else:
                        hdul[2].header["KS_P"] = (None, "p value of kstest vs noise")

                shutil.copy2(best_fit_gpath, tmp_fits_path_gname)
                
            except FileNotFoundError:
                pass
        
        # For when Simultaneous fitting fails we don't want to use that residual mask
        # for calculating the residual. I think everything else is handled
        #else:
            # TODO: MAY NOT NEED THIS ANY MORE FIXED IN RESIDUAL CALC
        #    shutil.copy2(pj(tmp_masks_dir, f"{gname}_star-rm.fits"), sf_masks_dir)
        #    replacement_sf_masks.append(pj(sf_masks_dir, f"{gname}_star-rm.fits"))
        
        # Alternatively, if not final_galfit_output.success and not galfit_output.success:
        if not exists(f"{tmp_fits_path_gname}"):
            print(f"{gname} completely failed!!!")
            # This file avoids conflating ones which haven't run and ones which have nothing to show for it
            touch_failed_fits(gname, tmp_fits_dir)
            # for debugging
            #print(final_out_text)
            print()
            #failed.append(gname)
            continue
        
        if sp(f"hostname").stdout.split(".")[0] == "bayonet-09":
                
            tmp_fits_obj = OutputFits(tmp_fits_path_gname, load_default = False)
            tmp_fits_obj.to_png(tmp_fits_path = tmp_fits_path_gname,
                                tmp_png_path  = tmp_png_path,
                                out_png_dir   = out_png_dir
                                #cleanup = False # TEMPORARY UNTIL IMAGEMAGICK WORKS AGAIN
                               )

            shutil.copy2(f"{pj(out_png_dir, gname)}_combined.png", pj(out_dir, gname))

        # No point in doing this in parallel because race conditions
        if not parallel and not use_async:
            for galfit_out in glob(pj(cwd, "galfit.*")):
                ext_num = galfit_out.split(".")[-1]
                shutil.move(galfit_out, pj(out_dir, gname, f"{gname}_galfit.{ext_num}"))
                
        # Residual calculation, now done all at the same time! And added to the FITS header
        if simultaneous_fitting:
            #_ = fill_objects(gname, 1, tmp_fits_dir, sf_masks_dir)
            
            # Remove copied masks
            rm_files(*replacement_sf_masks)

        shutil.copy2(tmp_fits_path_gname, pj(out_dir, gname, f"{gname}_galfit_out.fits"))
        
        # Creating a version of the galaxy for SpArcFiRe to use
        # Do this last just in case there's an issue
        _ = generate_model_for_sparcfire(
            tmp_fits_path_gname, 
            base_galfit_cmd, 
            in_dir
        )
        
        print()
        
    if aggressive_clean:
        #print("Aggressively cleaning... in, temp galfit output, and psf directories.")
        print("Aggressively cleaning... temp galfit output files.")
        # to_del_gal_in    = (pj(in_dir, f"{gname}.fits") for gname in galaxy_names
        #                     #if exists(pj(in_dir, f"{gname}.fits"))
        #                    )

        # Spell these out explicitly to avoid needing to glob anything... all that searching takes time!
        to_del = []
        for (b,d) in b_d_magnitudes:
        #for b in b_magnitudes:
            #prefix = f"m{b}m{d}"
            prefix = f"m{10*b:.0f}"
            
            # Parameter Search outputs
            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_galfit_out.fits") for gname in galaxy_names
            ])

            # Feedmes
            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}.in") for gname in galaxy_names
            ])
            
            # For if no arms
            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_bulge+disk.in") for gname in galaxy_names
            ])
            
            # Just doing all possibilities and letting rm_files
            # handle it if they don't exist
            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_bulge.in") for gname in galaxy_names
            ])

            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_bulge+disk_for_arms+arms.in") for gname in galaxy_names
            ])

            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_bulge+disk_for_arms+arms.in") for gname in galaxy_names
            ])

            to_del.extend([
                pj(tmp_fits_dir, f"{prefix}_{gname}_bulge+disk+disk_for_arms+arms.in") for gname in galaxy_names
            ])
            
        # Masks
        to_del.extend([
            pj(tmp_masks_dir, f"{gname}_star-rm.fits") for gname in galaxy_names
        ])
        
        # For sparcfire inputs
        to_del.extend([
            pj(tmp_fits_dir, f"{gname}_for_sparcfire.in") for gname in galaxy_names
        ])
        
        # Use try except instead of checking existence to save time
        # In most cases these files *should* exist
        rm_files(*to_del)
        print("Done!")
    
    return failed
        
if __name__ == "__main__":
    
    # cwd, in_dir, tmp_dir, out_dir, num_steps, rerun, verbose, *galaxy_names
    cmd_line_input = sys.argv[1:]
    
    kwargs_input = {}
    
    for kw in cmd_line_input:
        k, v = kw.split("=")
        v = v.strip()
        if v.lower() == "false":
            v = False
        elif v.lower() == "true":
            v = True
            
        kwargs_input[k] = v

    print(kwargs_input)
    _ = main(**kwargs_input)
    #write_failures(os.getcwd(), failures)
    
