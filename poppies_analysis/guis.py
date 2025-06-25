## Updated by Farhan Hasan (fhasan@stsci.edu)
## 6/5/25

#### This python script is mostly for DS9 GUIs

import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
import os
from pathlib import Path
from poppies_analysis import *
from distutils.sysconfig import *
from array import array
from glob import glob
import pandas as pd
from ast import literal_eval
import time


# POPPIES gui helpers
from guis_helpers import extract_image_extensions_key
from guis_helpers import display_image_in_DS9, display_images_in_DS9

# from samp_helpers import send_fits_to_ds9, send_regions_to_ds9 

from samp_helper import DS9SAMPHelper, get_ds9_samp_helper   ## FH added 6/2/25

# xpa resources:
# https://www.astro.louisville.edu/software/xmccd/archive/xmccd-4.1/xmccd-4.1e/docs/xpa/xpa.pdf
# http://ds9.si.edu/doc/ref/xpa.html#lock    

## Astropy SAMP resources:
# https://docs.astropy.org/en/stable/samp/index.html


# Global SAMP helper instance
_ds9_samp = None

def get_ds9_samp_helper():
    """Get or create DS9 SAMP helper instance."""
    global _ds9_samp
    if _ds9_samp is None or not _ds9_samp.connected:
        _ds9_samp = DS9SAMPHelper()
    return _ds9_samp

def cleanup_ds9_samp():
    """Clean up SAMP connection when done."""
    global _ds9_samp
    if _ds9_samp:
        _ds9_samp.disconnect()
        _ds9_samp = None

# Context manager for automatic cleanup
class DS9SAMPContext:
    """Context manager for DS9 SAMP operations."""
    
    def __enter__(self):
        return get_ds9_samp_helper()
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        cleanup_ds9_samp()


## FH updated 6/2/25 - converted to SAMP
def showSpec2D_POPPIES(parno, obid, filter = "F444W", path_to_data = "", zsc = "squared"):
    """Display spec2D cutouts in DS9 using SAMP.

    This can handle missed data in a given filter since there is a key
    to check against given the nature of how the spec2D files are setup.
    Frames will be made for missing data but will show nothing to indicate
    that the data is missing and that display layout is preserved.
    """
    ds9_title = "POPPIES_spec2D"

    try:
            
        # Get SAMP helper
        samp = get_ds9_samp_helper()

        #### Identify the correct client ID:
        try:             
            labels = samp.identify_client_labels()
            for i, (ids, labs) in enumerate(labels.items()):
                for j,name in enumerate(labs):
                    if name == ds9_title:
                        cl = j

            clientID = labels["client_id"][cl]
            
            # print('spec2d ', clientID)

        except Exception as e:
            print("Client ID not found due to error: ",e)


        def parse_filename(path_to_data, parno, obid):
            f1 = glob(path_to_data + f"/{parno}/jw*/A_{filter}_R_{obid}.V4.fits")  
            f2 = glob(path_to_data + f"/{parno}/jw*/A_{filter}_C_{obid}.V4.fits") 
            f3 = glob(path_to_data + f"/{parno}/jw*/B_{filter}_R_{obid}.V4.fits")  
            f4 = glob(path_to_data + f"/{parno}/jw*/B_{filter}_C_{obid}.V4.fits")         
            return f1,f2,f3,f4

        spec2D_files = parse_filename(path_to_data, parno, obid)
        
        # Check if given grism files exist or not:
        if len(spec2D_files[0]) > 0:
            grismr_file = spec2D_files[0][0]
        elif len(spec2D_files[2]) > 0: 
            grismr_file = spec2D_files[2][0]
        else:
            grismr_file = None

        if len(spec2D_files[1]) > 0:
            grismc_file = spec2D_files[1][0]
        elif len(spec2D_files[3]) > 0: 
            grismc_file = spec2D_files[3][0]
        else:
            grismc_file = None


        ### First for GRISM-R:
        if grismr_file is not None:

            fits_path = Path(grismr_file)

            # Extract the data information that we want to be displayed in ds9
            spec2D_key_DS9 = extract_image_extensions_key(grismr_file)

            # Display the results
            for key, _ in spec2D_key_DS9.items():
                # Get key-values for frame id and extension
                frame_id = int(spec2D_key_DS9[key]["frame_id"]) - 1    ## -1 bc first extension is 1D spec
                ext = spec2D_key_DS9[key]["ext"]

                # Display frame
                if ext is not None:

                    # Width and height of 2D spec trace
                    xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]
                    xpix, ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2.5*ypix_raw))

                    # print(ext,frame_id,grismr_file)

                    # Frame operations using SAMP
                    samp.send_command_to_ds9(f"frame {frame_id}", client_id=clientID)
                    samp.send_command_to_ds9("frame clear", client_id=clientID)

                    samp.send_command_to_ds9("frame refresh", client_id=clientID)
                    
                    # Load FITS file
                    samp.send_image_to_ds9(fits_path,grismr_file, client_id=clientID, extension=ext)
                    
                    # samp.send_command_to_ds9("frame center", client_id=clientID)
                    # samp.send_command_to_ds9("zoom 1", client_id=clientID)
        else:
            # Clear frames 1, 2, 3 if no GRISM-R file
            for frame_id in [1, 2, 3]:
                samp.send_command_to_ds9(f"frame {frame_id}", client_id=clientID)
                samp.send_command_to_ds9("frame clear", client_id=clientID)

        ### Now for GRISM-C:
        if grismc_file is not None:

            fits_path = Path(grismc_file)

            # Extract the data information that we want to be displayed in ds9
            spec2D_key_DS9 = extract_image_extensions_key(grismc_file)

            # Display the results
            for key, _ in spec2D_key_DS9.items():
                # Get key-values for frame id and extension
                frame_id = (int(spec2D_key_DS9[key]["frame_id"]) - 1) + 3    ## +3 bc new orientation here
                ext = spec2D_key_DS9[key]["ext"]

                # Display frame
                if ext is not None:
                    # Width and height of 2D spec trace
                    xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]
                    xpix, ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2.5*ypix_raw))

                    # Frame operations using SAMP
                    samp.send_command_to_ds9(f"frame {frame_id}", client_id=clientID)
                    samp.send_command_to_ds9("frame clear", client_id=clientID)

                    samp.send_command_to_ds9("frame refresh", client_id=clientID)
                    
                    # print(ext,frame_id,grismr_file)
                    
                    # Load FITS file
                    samp.send_image_to_ds9(fits_path,grismc_file, client_id=clientID, extension=ext)
                    
                    samp.send_command_to_ds9("frame center", client_id=clientID)
                    samp.send_command_to_ds9("zoom 1", client_id=clientID)
        else:
            # Clear frames 4, 5, 6 if no GRISM-C file
            for frame_id in [4, 5, 6]:
                samp.send_command_to_ds9(f"frame {frame_id}", client_id=clientID)
                samp.send_command_to_ds9("frame clear", client_id=clientID)


        # Format display into grid
        samp.send_command_to_ds9("tile", client_id=clientID)
        samp.send_command_to_ds9("tile grid mode manual", client_id=clientID)
        samp.send_command_to_ds9("tile mode row", client_id=clientID)
        samp.send_command_to_ds9("tile grid layout 1 6", client_id=clientID)

        # Set display dimensions
        try:
            samp.send_command_to_ds9(f"width {xpix}", client_id=clientID)
            samp.send_command_to_ds9(f"height {ypix}", client_id=clientID)
        except Exception as e:
            print('Could not correctly display spec2D, Reason: ', e)
            pass
        
        # Configure each frame
        # for fno in range(6, 0, -1):
        for fno in range(1,8):
            samp.send_command_to_ds9(f"frame {fno}", client_id=clientID)
            # samp.send_command_to_ds9(f"frame clear", client_id=clientID)
            samp.send_command_to_ds9("frame center", client_id=clientID)
            samp.send_command_to_ds9("zoom 1", client_id=clientID)
            
            # Display wavelength beside alpha
            samp.send_command_to_ds9("wcs fk5", client_id=clientID)
            samp.send_command_to_ds9("wcs skyformat degrees", client_id=clientID)
            samp.send_command_to_ds9(f"pan to {int(xpix_raw/2)} {int(ypix_raw/2)} image", client_id=clientID)

        # Set scaling and color map
        samp.send_command_to_ds9("scale mode zscale", client_id=clientID)
        samp.send_command_to_ds9("scale mode 99.5", client_id=clientID)
        samp.send_command_to_ds9(f"scale {zsc}", client_id=clientID)
        samp.send_command_to_ds9("cmap heat", client_id=clientID)

        # Lock settings
        samp.send_command_to_ds9("lock colorbar", client_id=clientID)
        samp.send_command_to_ds9("lock scale", client_id=clientID)
        # samp.send_command_to_ds9("zoom 1", client_id=clientID)
        samp.send_command_to_ds9(f"pan to 100 100 image", client_id=clientID)

    except Exception as e:
        print("Could not display spec2D for reason: ", e)
        return False


def find_file(directory, filename):
    """
    Scans a directory to find a file with the given name.

    Parameters:
    directory (str): The path of the directory to scan.
    filename (str): The name of the file to find.

    Returns:
    str: The path to the file if found, otherwise None.
    """
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None


## FH updated 5/30/25
def showDirect_POPPIES(parno, filter="F444W", path_to_data="", path_to_out=""):
    """Displays direct images for each of the filters"""

    ### KVN's quick fix to images having different names in different fields
    # specify the images to be displayed in DS9
    # grism_file = glob(path_to_data + '/POPPIES'+str(parno)+'/DATA/*gr150*_drz_sci.fits')
    ### FH updated 1/6/25
    # grism_file = glob(path_to_data + '/POPPIES'+str(parno)+'/DATA/*A_F444W_i2d.fits')
    # grism_file = str(grism_file[0]).split(path_to_data)[1]
    # grism_file_2 = str(grism_file).split('DATA/')[1]

    image_files = glob(path_to_data + str(parno)+'/*_{}_i2d.fits'.format(filter))
    
    image_file = image_files[0]

    # image_file_ext = str(image_file).split('DirectImages/')
    image_file_ext = str(image_file).split(str(parno)+'/')
    image_file_2 = image_file_ext[1]

    images = {
        f"{filter}": [
            f"{image_file_2}"
            # f"POPPIES{parno}_{grism_file_ext}_f115w-gr150c_drz_sci.fits",
            # f"POPPIES{parno}_{grism_file_ext}_f115w-gr150r_drz_sci.fits",
        ]
    }

    #### KVN assuming some path structure because this doesn't work for me...
    ## FH updated 1/6/25
    # find the full paths of the direct images
    image_paths = {}
    for filter_name, filenames in images.items():
        paths = []
        for filename in filenames:
            full_path = find_file(path_to_data, filename)
            # print(full_path)
            if full_path:
                paths.append(full_path)
            else:
                print(f"Warning: File {filename} not found; assuming it lives in path_to_data/#/")
                full_path = find_file(path_to_data+str(parno)+'/', filename)
                paths.append(full_path)
                # print(paths)
                
        image_paths[filter_name] = paths

    
    # specify the region files to show along with images, if any.
    # region_files = {
    #     "f115w": [f"POPPIES{parno}regions_phot.reg", f"POPPIES{parno}F115c_grism.reg", f"POPPIES{parno}F115r_grism.reg"],
    #     "f150w": [f"POPPIES{parno}regions_phot.reg", f"POPPIES{parno}F150c_grism.reg", f"POPPIES{parno}F150r_grism.reg"],
    #     "f200w": [f"POPPIES{parno}regions_phot.reg", f"POPPIES{parno}F200c_grism.reg", f"POPPIES{parno}F200r_grism.reg"],
    # }
    region_files = {
        f"{filter}": [f"POPPIES{parno}regions_phot_{filter}.reg"]
    }

    # find the full paths of the region files
    region_paths = {}
    for filter_name, filenames in region_files.items():
        paths = []
        for filename in filenames:
            full_path = find_file(path_to_out, filename)
            # print(full_path)
            if full_path:
                paths.append(full_path)
            else:
                # KVN: assuming same pathing structure
                print(f"Warning: File {filename} not found; assuming it lives in path_to_out/#/")
                full_path = find_file(path_to_out + str(parno)+ '/', filename)
                paths.append(full_path)
                
        region_paths[filter_name] = paths

    # display_images_in_DS9(image_paths, region_files=region_paths)
    display_images_in_DS9(image_paths, regions=region_paths, tile_mode=False)

    return


## FH 6/2/25 
def panDirect_POPPIES(ra, dec, zoom=4):
    # print(ra, dec)
    
    ds9_title = "POPPIES_DIRECT"

    samp = get_ds9_samp_helper()

    ### Identify the correct client ID:
    try:             
        labels = samp.identify_client_labels()
        for i, (ids, labs) in enumerate(labels.items()):
            # print(i,ids,labs)
            for j,name in enumerate(labs):
                if name == ds9_title:
                    cl = j

        clientID = labels["client_id"][cl]
        # print(clientID)

        # samp.send_command_to_ds9("frame refresh", client_id=clientID)
        
        # samp.send_command_to_ds9("frame 1", client_id=clientID)

        samp.send_command_to_ds9(f"pan to {ra[0]} {dec[0]} wcs degrees", client_id=clientID)

        # samp.send_command_to_ds9("pan to 1000000 10000000 image", client_id=clientID)

        samp.send_command_to_ds9(f"zoom to {zoom}", client_id=clientID)
        # samp.send_command_to_ds9("frame lock wcs", client_id=clientID)

        # samp.send_command_to_ds9("frame 1", client_id=clientID)
        # samp.send_command_to_ds9(f"frame clear", client_id=clientID)
        # samp.send_command_to_ds9("frame center", client_id=clientID)
        # samp.send_command_to_ds9("zoom 1", client_id=clientID)
        
        # Display wavelength beside alpha
        # samp.send_command_to_ds9("wcs fk5", client_id=clientID)
        # samp.send_command_to_ds9("wcs skyformat degrees", client_id=clientID)
        # samp.send_command_to_ds9("pan to 100 100 image", client_id=clientID)


    except Exception as e:
        print("Client ID not found due to error: ",e)


    return



#### XPA versions:


# ## FH updated 1/24/25
# def panDirect_POPPIES(ra, dec):
#     # print(ra, dec)

#     ds9_title = "POPPIES_DIRECT"
#     # for fno in [1,4,7]:

#     ## will do loop when more filters are added
#     # for fno in [1]:
#     cmd = f"xpaset -p {ds9_title} frame 1"  # + str(fno)
#     os.system(cmd)
#     cmd = f"xpaset -p {ds9_title} pan to {ra[0]} {dec[0]} wcs degrees"
#     os.system(cmd)

#     # # zoom
#     # cmd = f"xpaset -p {ds9_title} zoom 2"
#     # os.system(cmd)

#     # os.system(f"xpaset -p {ds9_title} zoom in")
#     os.system(f"xpaset -p {ds9_title} zoom to 4")


# ## FH updated 5/30/25
# def showSpec2D_POPPIES(parno, obid, filter = "F444W", path_to_data = "", zsc = "squared"):
#     """Display spec2D cutouts in DS9.

#     This can handle missed data in a given filter since there is a key
#     to check against given the nature of how the spec2D files are setup.
#     Frames will be made for missing data but will show nothing to indicate
#     that the data is missing and that display layout is preserved.
#     """

#     # try:
#     #     secats = glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/DIRECT_GRISM/POPPIES*phot*.fits")  # MDR 2022/05/17
#     # except:
#     #     secats = glob(path_to_data + "/POPPIES" + str(parno) + "/Products/POPPIES*phot*.fits")  # KVN allowing for different path structure (?)


#     def parse_filename(path_to_data, parno, obid):
#         # return path_to_data + f"/POPPIES{parno}/spec2D/POPPIES{parno}_{obid:05d}.2D.fits"
#         f1 = glob(path_to_data + f"/{parno}/jw*/A_{filter}_R_{obid}.V4.fits")  
#         f2 = glob(path_to_data + f"/{parno}/jw*/A_{filter}_C_{obid}.V4.fits") 

#         f3 = glob(path_to_data + f"/{parno}/jw*/B_{filter}_R_{obid}.V4.fits")  
#         f4 = glob(path_to_data + f"/{parno}/jw*/B_{filter}_C_{obid}.V4.fits")         
        
#         return f1,f2,f3,f4


#     spec2D_files = parse_filename(path_to_data, parno, obid)

#     #check if given grism 1d files exist or not:
#     if len(spec2D_files[0]) > 0:
#         grismr_file = spec2D_files[0][0]

#     elif len(spec2D_files[2]) > 0: 
#         grismr_file = spec2D_files[2][0]

#     else:
#         grismr_file = None

#     if len(spec2D_files[1]) > 0:
#         grismc_file = spec2D_files[1][0]

#     elif len(spec2D_files[3]) > 0: 
#         grismc_file = spec2D_files[3][0]

#     else:
#         grismc_file = None

#     ### first for GRISM-R:
#     # spec2D_file = parse_filename(path_to_data, parno, obid)[0][0]

#     SPEC2D_TITLE_DS9 = "POPPIES_spec2D"


#     if grismr_file is not None:
            
#         # for a given data file extract the data information that we want to be displayed in ds9
#         spec2D_key_DS9 = extract_image_extensions_key(grismr_file)


#         # display in ds9 instance of a given "title"

#         # now display the results.
#         for key, _ in spec2D_key_DS9.items():

#             # get key-values for frame id and extension
#             frame_id = int(spec2D_key_DS9[key]["frame_id"]) - 1    ## -1 bc first extension is 1D spec
#             ext = spec2D_key_DS9[key]["ext"]

#             # # display frame
#             if ext is not None:
                
#                 #width and height of 2D spec trace
#                 xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]

#                 xpix,ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2*ypix_raw))

#                 # frame id
#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}"

#                 os.system(command)

#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

#                 os.system(command)

#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} fits {grismr_file}[{ext}]"

#                 os.system(command)

#                 os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
#                 os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")
                
#         # else:
#         #     command = f"xpaset -p {SPEC2D_TITLE_DS9} frame clear"
#         #     os.system(command)
    
#     else:
#         for frame_id in [1,2,3]:
#             # frame id
#             os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}")
#             os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame clear")

#     # # formating display into grid
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile")
#     # # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid mode manual")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile mode row")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid layout 1 6")

#     # # config display properties
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 1")
#     # os.system("xpaset -p {SPEC2D_TITLE_DS9} lock frame physical")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock colorbar")

#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode user")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale yes")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits yes")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap grey")


#     # for fno in [1,2,3]:
#     #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " +str(fno))
#     #     # OLD: 39.5 and 32.5 are just 1/2 cutout length & 1/2 cutout width. Will need to be changed for NIRCam
#     #     # The length is set to the length of the cutout -8 (for 4 pix on each end)
#     #     # For NIRCam F444W, image size is 1408 x 200
#     #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 710 105 100 10# color=green} ")        
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale squared")
#     #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale")
#     #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits yes")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
#         # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1.5")

#         # command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

#         # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh")

#     ### Now for GRISM-C:
#     # spec2D_file = parse_filename(path_to_data, parno, obid)[1][0]


#     if grismc_file is not None:
#         # for a given data file extract the data information that we want to be displayed in ds9
#         spec2D_key_DS9 = extract_image_extensions_key(grismc_file)

#         # now display the results.
#         for key, _ in spec2D_key_DS9.items():

#             # get key-values for frame id and extension
#             frame_id = (int(spec2D_key_DS9[key]["frame_id"]) - 1) + 3    ## +3 bc new oreintation here
#             ext = spec2D_key_DS9[key]["ext"]

#             # xpix_vals,ypix_vals = [], []

#             # # display frame
#             if ext is not None:
                
#                 # print('GC ', frame_id)

#                 #width and height of 2D spec trace
#                 xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]

#                 xpix,ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2*ypix_raw))

#                 # frame id
#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}"

#                 os.system(command)

#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

#                 os.system(command)

#                 command = f"xpaset -p {SPEC2D_TITLE_DS9} fits {grismc_file}[{ext}]"

#                 os.system(command)

#                 os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
#                 os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")
#                 # xpix_vals.append(xpix),ypix_vals.append(ypix)


#             # else:
#             #     command = f"xpaset -p {SPEC2D_TITLE_DS9} frame clear"
#             #     os.system(command)


#     else:
#         for frame_id in [4,5,6]:
#             # frame id
#             os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}")
#             os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame clear")

#     # formating display into grid
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid mode manual")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile mode row")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid layout 1 6")

# #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width 1600")
# #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height 1200")


#     try:  #FH modified 3/6/25
#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")
#     except Exception as e:
#         print('Could not correctly display spec2D, Reason: ', e)
#         pass

#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale {zsc}")

#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock colorbar")
#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")

#     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")

#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 7")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")

#     # # for fno in np.arange(6,0,-1):

#     for fno in np.arange(1,8,1):

#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " + str(fno))
#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")

#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")
        
#         # FH 2/17/25: display wavelength beside alpha
#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} wcs fk5")

#         os.system(f"xpaset -p {SPEC2D_TITLE_DS9} wcs skyformat degrees")


#     #     # For NIRCam F444W, image size is 1408 x 200
#     #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 710 105 100 10# color=green} ")        
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")
#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale {zsc}")

#     #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

#         # ## Need a 180 degree rotation for R-spectra:
#         # if fno in [1,2,3]:
#         #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} rotate 180")

#         # if fno==7:
#         #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")



#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 7")
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")

#     # # Go to the first frame:    
#     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 1")
       