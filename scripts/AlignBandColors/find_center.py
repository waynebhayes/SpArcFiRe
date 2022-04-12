import os
import subprocess
import time

import numpy as np
from astropy.modeling import functional_models

from galaxy import Star


class CurveFitError(Exception):
    def __init__(self, message):
        self.message = message


def estimate_center(img : "ndarray", star : Star, percent_img_to_explore : float = 0.025) -> Star:
    """
    Given an image and an approximate star location, will fit a surface to the star and output the results

    Arguments:
        img                       (ndarray) : Image containing the star
        star                      (Star)    : Star object containing the approximate star position
        percent_of_img_to_explore (float)   : Percent of image to explore to look for the star (since SourceExtractor's
                                              points are often not centered, we have to center the image ourself)
    
    Returns:
        Star : Star object with the updated gamma, alpha, x, and y values found from the fit
    """
    s_size = 7
    point = (int(star.x), int(star.y))
    h, w = img.shape
    img_dist = int(h * percent_img_to_explore)

    # Calculate the x-y points for the zoomed up image
    y_min = max(0, point[1] - img_dist)
    y_max = min(h, point[1] + img_dist)
    x_min = max(0, point[0] - img_dist)
    x_max = min(w, point[0] + img_dist)
    
    sub_img = img[y_min : y_max, x_min : x_max]

    # Brightness pixel of the zoomed up image (roughly the center of the star)
    max_pt = np.unravel_index(np.argmax(sub_img), sub_img.shape)
    max_pt = (max_pt[1] + x_min, max_pt[0] + y_min)
 
    dist = np.sqrt((point[0] - max_pt[0]) ** 2 + (point[1] - max_pt[1]) ** 2)
    
    if dist > 10:
        raise CurveFitError(f"Distance between sextractor point and max point is too big: {dist} pixels.")

    # Zoom up on the image so that the brightest pixel is in the center (should give a better fit)
    ys_min = max(0, max_pt[1] - s_size)
    ys_max = min(h, max_pt[1] + s_size)
    xs_min = max(0, max_pt[0] - s_size)
    xs_max = min(w, max_pt[0] + s_size)
    
    centered_img = img[ys_min : ys_max, xs_min : xs_max]
    
    x = "["
    for i in range(centered_img.shape[1]):
        for j in range(centered_img.shape[0]):
            x += f"[{j}, {i}], "
    x = x[:-2] + "]"

    y = "["
    for value in centered_img.ravel():
        y += f"{value}, "
    y = y[:-2] + "]"
     
    max_pt = np.unravel_index(np.argmax(centered_img), centered_img.shape)
 
    # Input parameter Amplitude, x0, y0, gamma, alpha
    params = '[{}, {}, {}, 1.0, 1.0]'.format(str(centered_img[max_pt]), str(max_pt[1]), str(max_pt[0]))
    lb = '[0.0, {}, {}, 0.0, 0.0]'.format(str(max_pt[1] - np.sqrt(2)), str(max_pt[0] - np.sqrt(2)))
    ub = '[1000.0, {}, {}, 100.0, 100.0]'.format(str(max_pt[1] + np.sqrt(2)), str(max_pt[0] + np.sqrt(2)))
    
    abc_path   = os.path.abspath(os.path.dirname(__file__))
    curve_path = os.path.join(abc_path, "CurveFit/runFit")
    proc = subprocess.Popen([curve_path, x, y, params, lb, ub], stdout = subprocess.PIPE)
    
    timeout = 20
    while proc.poll() is None:
        time.sleep(1)
        timeout -= 1
        if timeout == 0:
            proc.terminate()
            raise CurveFitError(f"Ran out of time fiting star at {point}")
        
    out, err = proc.communicate()
    
    if proc.wait() != 0:
        raise CurveFitError(f"Error fitting curve to star at {point}")
   
    p = np.array(out.decode("utf-8").split(' ')).astype(float)
    f = functional_models.Moffat2D(amplitude = p[0], x_0 = p[1] + xs_min, y_0 = p[2] + ys_min, gamma = p[3], alpha = p[4])
    return Star(f.x_0.value, f.y_0.value, gamma = f.gamma.value, alpha = f.alpha.value, class_prob = star.class_prob)
