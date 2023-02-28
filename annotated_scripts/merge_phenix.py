import os  # provides functions for interacting with the operating system
import re  # provides regular expression matching operations
import sys  # provides various functions and variables that are used to
            # manipulate different parts of the Python runtime environment
import numpy as np  # multidimensional array objects and a collection of
                    # routines for processing of array
from PIL import Image  # Python Imaging Library; Image module provides a number
                       # of factory functions, including functions to load
                       # images from files, and to create new images.

# If number of inputs isn't equal to 3, print
# "Usage: <script_name> <input_folder> <destination_folder>" and exit with an
# error
if len(sys.argv) != 3:
    print("""Usage: %s <input_folder> <destination_folder>""" % sys.argv[0])
    exit(1)

# _ is the script directory, infolder is input folder, and outfolder is output
# folder. _ is used because it's a throwaway variable. map(function, iterable)
# extracts each of these paths.
# os.path.abspath() returns a normalized absolutized version of the pathname
# path
_, infolder, outfolder = map(os.path.abspath, sys.argv)
#phenix_re = re.compile('^r(\d{2})c(\d{2})f(\d{2})p(\d{2})-ch(\d)sk*.tiff$')

# Check if the input folder exists, and exit the program if not
if not os.path.exists(infolder):
    print("Input folder does not exist")
    exit(1)

# Check if the output folder exists, and create it if not
if not os.path.exists(outfolder):
    print('Creating outfolder', outfolder)
    os.makedirs(outfolder)

save_ch3 = False

# os.walk() generates filenames in a directory tree by walking it top-down (or
# bottom-up). For each directory it yields a three-item tuple
# (directory_[path, folder_names], [file_names]) or (path, [folder], [files]) in this case.
for path, folders, files in os.walk(infolder):
    # sorted() sorts from smallest to largest
    # lambda is used to filter for .tiff or .tif files and don't start with '.'
    # and contain 'ch1' in their name in files
    tiffFiles = sorted(filter(lambda x: ((x.endswith('.tiff') or x.endswith('.tif')) and not x.startswith('.') and 'ch1' in x), files))

    # if tiffFiles is empty, continue
    if not tiffFiles:
        continue

    # replace method used on the string path; input folder path replaced with
    # output folder path
    output_path = path.replace(infolder, outfolder)
    # if output path doesn't exist, make it (didn't we already do this??)
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)
    print(output_path)

    all_tiff = {}
    # For each .tiff file, create its full path by combining the .tiff file name
    # to the path string, then use that to open ch1 or ch2 images. Put ch2 in
    # its own list (why?)
    for tiff in tiffFiles:
        tiffpath = os.path.join(path, tiff)
        image_ch1 = Image.open(tiffpath)
        image_ch2 = Image.open(tiffpath.replace('ch1', 'ch2'))
        append_images = [image_ch2]

        # I'm assuming in some cases there's a third channel, and you'd have to
        # indicate this when submitting jobs like save_ch3=TRUE (how?) Either
        # way, open ch3 image and put it in the same list as ch2
        if save_ch3:
            path_ch3 = tiffpath.replace('ch1', 'ch3')
            image_ch3 = Image.open(path_ch3)
            append_images.append(image_ch3)
        # create new filename for merged image by joining output path with first
        # 12 characters of .tiff file (eliminates channel information)
        merged_filename = os.path.join(output_path, "%s.tiff" % tiff[:12])
        # use save() method to combine ch1 image with everything in the ch2/ch3
        # list and save it under merged_filename
        # save_all: save all frames of image instead of just the first one
        # compression: string specifying compression method (if no compression,
        # None)
        # append_images: list of images to append as additional frames
        # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html
        image_ch1.save(merged_filename, save_all=True, compression="None", append_images=append_images)
        print('Saved %s' % merged_filename)
print('DONE')

