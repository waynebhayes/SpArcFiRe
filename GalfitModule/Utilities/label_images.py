# Credit for underlying functionality goes to:
# https://www.blog.pythonlibrary.org/2021/02/16/creating-an-image-viewer-with-pysimplegui/

# USAGE: python label_images.py BASE_DIR DIRNAME_FOR_MISLABELED
# BASE_DIR is the directory that houses the pre-labeled images
# in our case, these are the images that pass/fail our verification tests
#
# DIRNAME_FOR_MISLABELED is the directory **inside** BASE_DIR which
# will hold the images which were mislabeled by the verification tests
# i.e. false positive/false negative
#
# On running, an interaction window will pop up with the buttons,
# prev, next, mislabeled. These allow the user to view images in 
# BASE_DIR (jpg only) and either move them to DIRNAME_FOR_MISLABLED 
# via 'mislabled' or continue viewing images via next/prev.
# hint: mogrify -format jpg *.png
# 
# Once all images have been viewed/moved, the window will close.
#
# Known bug: no images on initial run. Must click 'next' to see the second
# image in the set, then 'prev' to view the first before proceeding onwards.
# This is not a high priority but will be fixed eventually.

import glob
import PySimpleGUI as sg
from shutil import move
from os.path import exists
from os.path import join as pj
import sys

from PIL import Image, ImageTk


def parse_folder(path):
    images = glob.glob(f'{path}/*.jpg')# + glob.glob(f'{path}/*.png')
    return images

def load_image(path, window):
    try:
        image = Image.open(path)

        # TODO: Detect if horizontal or vertical and choose accordingly
        # Horizontal
        #image.thumbnail((800, 375))

        # Vertical
        image.thumbnail((375, 800))

        photo_img = ImageTk.PhotoImage(image)
        window["image"].update(data=photo_img)
    except:
        print(f"Unable to open {path}!")
        

def main():
    
    elements = [
        [sg.Image(key="image")],
        #[
            #sg.Text("Image File"),
            #sg.Input(size=(25, 1), enable_events=True, key="file"),
            #sg.FolderBrowse(),
        #],
        [
            sg.Button("Prev"),
            sg.Button("Next"),
            sg.Button("Mislabeled")
            #sg.Button("Non-Spiral")
        ]
    ]

    # Horizontal
    #window = sg.Window("Image Viewer", elements, size=(900, 400))
    
    # Vertical
    window = sg.Window("Image Viewer", elements, size=(400, 900))
    
    images = []
    location = 0

    img_path        = sys.argv[1]
    failure_path    = sys.argv[2] #pj(values["file"], "false_positive")
#    non_spiral_path = pj(basepath, sys.argv[3]) #pj(values["file"], "false_positive")
    images = parse_folder(img_path)#values["file"])
    while images:
        load_image(images[location], window)
        event, _ = window.read()

        if event == "Exit" or event == sg.WIN_CLOSED:
            break

        if event == "Next" and images:
            if location == len(images) - 1:
                break
            else:
                location += 1
            continue

            #load_image(images[location], window)

        if event == "Prev" and images:
            if location == 0:
                #location = len(images) - 1
                pass
            else:
                location -= 1
            #load_image(images[location], window)
            continue

        if event == "Mislabeled" and images:
            if exists(images[location]):
                print(f"Moving {images[location]} to {failure_path}")
                move(images[location], failure_path)
                images.pop(location)
            else:
                print("image has already been moved... continuing")

            if location == len(images) - 1:
                break
            #else:
            #    location += 1

            #load_image(images[location], window)

#        if event == "Non-Spiral" and images:
#            if exists(images[location]):
#                print(f"Moving {images[location]} to {non_spiral_path}")
#                move(images[location], non_spiral_path)
#                images.pop(location)
#            else:
#                print("image has already been moved... continuing")

#            if location == len(images) - 1:
#                break
#            else:
#                location += 1

#            load_image(images[location], window)

    window.close()


if __name__ == "__main__":
    main()