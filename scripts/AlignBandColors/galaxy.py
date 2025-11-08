from copy import copy
from collections import OrderedDict

import numpy as np

import load_gals


# Exception is raised when a waveband color is requested that is not present
class InvalidGalColorError(Exception): pass

# Exception is raised when an error is encountered when cropping the images
class CroppingError(Exception): pass

# Exception is raised when creating a galaxy and wavebands do not match
class GalsAndStarsDoNotContainTheSameWavebandsError(Exception): pass

# Exception is raised when image dimensions are not the same across all wavebands loaded
class ImageDimensionsInconsistentError(Exception): pass


class Galaxy:
    """ Class that stores all information about a galaxy in all wavebands """
    def __init__(self, gal_dict, stars_dict, star_class_perc, name):
        # It is expected that gal_dict and stars_dict contain the same
        #    wavebands in the same order.

        if gal_dict.keys() != stars_dict.keys():
            raise GalsAndStarsDoNotContainTheSameWavebandsError
        
        self.gal_dict = gal_dict
        self.all_stars_dict = copy(stars_dict)
        
        # Keep a separate list of those that are likely stars
        self.stars_dict = OrderedDict()
        for color in stars_dict.keys():
            self.stars_dict.update({color: [s for s in stars_dict[color] if s.class_prob >= star_class_perc]})

        self.name = name

        # Check all dimensions are the same
        widths  = {gal[0].data.shape[1] for gal in self.gal_dict.values()}
        heights = {gal[0].data.shape[0] for gal in self.gal_dict.values()}
        if len(widths) != 1 or len(heights) != 1:
            raise ImageDimensionsInconsistentError
        
        self.width  = widths.pop()
        self.height = heights.pop()
        self.num_wb = len(self.gal_dict)
        
    
    def images(self, color : str = None) -> "ndarray":
        """ 
        Gets the galaxy image associated with the given color

        Arguments:
            color (str) : Value in griuz
        
        Returns:
            Tuple of all galaxy color - image pairs if color is None, otherwise just the given color image
        """
        if color is None:
            return [(color, img[0].data) for color, img in self.gal_dict.items()]

        elif color in self.gal_dict.keys():
            return self.gal_dict[color][0].data
        
        raise InvalidGalColorError

    
    def crop_images_to_galaxy(self) -> None:
        """
        Runs source extractor on the images and crops the images
        down to only the galaxy (modifies the class images!)
        """

        left, right, top, bottom = 0, 0, 0, 0
        
        for _, img in self.images():
            seg_img = load_gals.get_seg_img(img)
            gal_val = seg_img[int(seg_img.shape[0] / 2), int(seg_img.shape[1] / 2)]
            inds = np.argwhere(seg_img == gal_val)
            x_inds, y_inds = inds[:,1], inds[:,0]
            left += np.min(x_inds); right += np.max(x_inds)
            top += np.min(y_inds); bottom += np.max(y_inds)
        
        left /= float(self.num_wb); right /= float(self.num_wb)
        top /= float(self.num_wb); bottom /= float(self.num_wb)

        center_x, center_y = self.width / 2, self.height / 2
        size = int(max(center_x - left, right - center_x, center_y - top, bottom - center_y)) + 2
        left, right = int(center_x - size), int(center_x + size)
        top, bottom = int(center_y - size), int(center_y + size)

        # Make sure the values found are valid
        try:
            assert top >= 0; assert left >= 0; 
            assert bottom <= self.height; assert right <= self.width
        except AssertionError:
            size -= 2
            left, right = int(center_x - size), int(center_x + size)
            top, bottom = int(center_y - size), int(center_y + size)
            if top < 0 or left < 0 or bottom > self.height or right > self.width:
                raise CroppingError

        # Crop the images
        for c in self.gal_dict.keys():
            self.gal_dict[c][0].data = self.gal_dict[c][0].data[top:bottom, left:right]
        
        return left, right, top, bottom


    def stars(self, color : str = None) -> "[Star]": 
        """ 
        Gets the stars associated with the given color

        Arguments:
            color (str) : Value in griuz
        
        Returns:
            Tuple of all color - star pairs if color is None, otherwise just the given color's stars
        """
        if color is None:
            return self.stars_dict.values()
        
        elif color in self.stars_dict.keys():
            return self.stars_dict[color]

        raise InvalidGalColorError


    def gen_img_star_pairs(self) -> "(color, image, stars)":
        """
        Generator to retrieve color-galaxy-star tuples
        """
        for color, gal, stars in zip(self.gal_dict.keys(), self.gal_dict.values(), self.stars_dict.values()):
            yield (color, gal[0].data, stars)   


    def colors(self) -> "[Colors]":
        """
        Returns all colors present
        """
        return self.gal_dict.keys()

    
    def add_borders(self, b_size : int) -> None:
        """
        Adds a black border of b_size pixels around all images
        
        Arguments:
            b_size (int) : size of border in pixels
        """
        for img in self.gal_dict.values():
            img[0].data = np.pad(img[0].data, b_size, 'constant')
            img[0].header['CRPIX1'] += b_size
            img[0].header['CRPIX2'] += b_size
    
    def close(self) -> None:
        """
        Closes all of the astropy.fits objects
        """
        for gal in self.gal_dict.values():
            gal.close()



class Star:
    """ Stores all of the information about a star from Source Extractor and its fit"""    
    def __init__(self, x, y, class_prob = None, gamma = None, alpha = None):
        self.x = x
        self.y = y
        self.class_prob = class_prob
        self.gamma = gamma
        self.alpha = alpha

    def info(self):
        return f"({self.x}, {self.y}, {self.gamma}, {self.alpha}, {self.class_prob})"
    
    def __str__(self):
        return f"({self.x}, {self.y})"

    def __sub__(self, star) -> "ndarray":
        """ Overload for subtracting stars (does axis-wise subtraction) """
        return np.array((self.x - star.x, self.y - star.y))
