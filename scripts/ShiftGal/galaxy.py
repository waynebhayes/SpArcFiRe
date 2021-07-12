from collections import OrderedDict
import numpy as np
from copy import copy
import load_gals


class InvalidGalColorError(Exception): pass

class CroppingError(Exception): pass

class GalsAndStarsDoNotContainTheSameWavebandsError(Exception): pass


class Galaxy:

    def __init__(self, gal_dict, stars_dict, star_class_perc, name):
        # It is exptected that gal_dict and stars_dict contain the same
        # wavebands in the same order.

        if gal_dict.keys() != stars_dict.keys():
            raise GalsAndStarsDoNotContainTheSameWavebandsError
        
        self.gal_dict = gal_dict
        self.all_stars_dict = copy(stars_dict)
        # keep a separate list of those that are likely stars
        self.stars_dict = OrderedDict()
        for color in stars_dict.keys():
            self.stars_dict.update({color: [s for s in stars_dict[color] if s.class_prob >= star_class_perc]})

        self.name = name
        if 'g' in self.gal_dict.keys():
            self.width = self.gal_dict['g'][0].data.shape[1]
            self.height = self.gal_dict['g'][0].data.shape[0]
        self.num_wb = len(self.gal_dict)
        
    
    def images(self, color = None):
        if color is None:
            return [(color, img[0].data) for color, img in self.gal_dict.iteritems()]

        elif color in self.gal_dict.keys():
            return self.gal_dict[color][0].data
        
        raise InvalidGalColorError

    
    def crop_images_to_galaxy(self):
        """Runs source extractor on the images and crops the images
           down to only the galaxy"""

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

        # make sure the values found are valid
        try:
            assert top >= 0; assert left >= 0; 
            assert bottom <= self.height; assert right <= self.width
        except AssertionError:
            size -= 2
            left, right = int(center_x - size), int(center_x + size)
            top, bottom = int(center_y - size), int(center_y + size)
            if top < 0 or left < 0 or bottom > self.height or right > self.width:
                raise CroppingError

        # crop the images
        for c in self.gal_dict.keys():
            self.gal_dict[c][0].data = self.gal_dict[c][0].data[top:bottom, left:right]
        
        return left, right, top, bottom


    def stars(self, color = None): 
        if color is None:
            return self.stars_dict.values()
        
        elif color in self.stars_dict.keys():
            return self.stars_dict[color]

        raise InvalidGalColorError


    def gen_img_star_pairs(self):
        for color, gal, stars in zip(self.gal_dict.keys(), self.gal_dict.values(), self.stars_dict.values()):
            yield (color, gal[0].data, stars)   


    def colors(self):
        return self.gal_dict.keys()

    
    def add_borders(self, b_size):
        for img in self.gal_dict.values():
            img[0].data = np.pad(img[0].data, b_size, 'constant')
            img[0].header['CRPIX1'] += b_size
            img[0].header['CRPIX2'] += b_size
    
    def close(self):
        for gal in self.gal_dict.values():
            gal.close()


class Star:

    def __init__(self, x, y, class_prob = None, gamma = None, alpha = None):
        self.x = x
        self.y = y
        self.class_prob = class_prob
        self.gamma = gamma
        self.alpha = alpha

    def info(self):
        return '({}, {}, {}, {}, {})'.format(self.x, self.y, self.gamma, self.alpha, self.class_prob)
    
    def __str__(self):
        return '({}, {})'.format(self.x, self.y)

    def __sub__(self, star):
        return np.array((self.x - star.x, self.y - star.y))
