## Updated by Farhan Hasan (fhasan@stsci.edu)
## 6/5/25

#!/usr/bin/env python3
"""
Updated guis_helpers.py functions using SAMP instead of XPA
Migration from XPA-based DS9 communication to SAMP-based communication
"""

#### This python script is mostly for DS9 GUIs (as a sort of secondary script)

import os
import time
import numpy as np
from pathlib import Path
from astropy.samp import SAMPIntegratedClient
from astropy.io import fits
import logging
from samp_helper import DS9SAMPHelper, get_ds9_samp_helper   ## FH added 6/2/25


# Set up logging
logger = logging.getLogger(__name__)

# logging.basicConfig(level=logging.INFO,filemode='w')

# logfile = str('/Users/fhasan/Desktop/Research_STScI/POPPIES/samptests.log')

# # Set log file output
# handler = logging.FileHandler(logfile,'w')
# handler.setLevel(logging.INFO)

# # Create a logging format
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)

# # Add the handlers to the logger
# logger.addHandler(handler)


### SAMP Helper stuff:

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


### FH updated 1/6/25
def extract_image_extensions_key(filename):

    # frame number ids
    frame_ids = np.arange(2, 7)

    filter_orientations = [
        # ("F115W,199.0", "SCI"),
        # ("F115W,199.0", "CONTAM"),
        # ("F115W,290.0", "SCI"),
        # ("F115W,290.0", "CONTAM"),
        # ("F115W", "COMB"),
        # ("F150W,199.0", "SCI"),
        # ("F150W,199.0", "CONTAM"),
        # ("F150W,290.0", "SCI"),
        # ("F150W,290.0", "CONTAM"),
        # ("F150W", "COMB"),
        # ("F200W,199.0", "SCI"),
        # ("F200W,199.0", "CONTAM"),
        # ("F200W,290.0", "SCI"),
        # ("F200W,290.0", "CONTAM"),
        # ("F200W", "COMB"),
        # ("2D_COUNT"),
        # ("2D_CCOUNT"),
        # ("2D_MCOUNT"),
        # ("2D_ERR"),
        # ("2D_TCOUNT"),
        # ("2D_WAVE")
        ("2D_COUNT"),
        ("2D_CCOUNT"),
        ("2D_MCOUNT"),
        ("2D_ERR"),
        ("2D_TCOUNT"),
        ("2D_WAVE")        
        ]
    
    # filter + orientation of COUNT, CCOUNT, and MCOUNT images and their associated
    # frame number use to plot in DS9.
    spec2D_key = {
        key: {"frame_id": val, "ext": None, "x_pix": None, "y_pix": None}
        for key, val in zip(filter_orientations, frame_ids)
    }

    with fits.open(filename) as hdu:
        n_ext = len(hdu)

        extlist = [2,5,6]

        # loop over extensions to search for to identify which extension belongs to
        # which frame id in spec2D_key
        
        ext_out = 0 # running counter for frame no.

        for ext in range(2, n_ext-1):

            # extract header information - if one of the extensions we care about
            if ext in extlist:
                ext_out += 1
                header = hdu[ext].header

                # check if name is in desired reference names
                # if header["EXTNAME"] in ["2D_COUNT"]:

                extname = header["EXTNAME"]

                naxis1, naxis2 = header["NAXIS1"], header["NAXIS2"]

                # extver = header["EXTVER"].split(",")
                # extver = "F444W"

                # if extver is on the filter name w/out orientation
                # its the combined sci images images so label COMB
                # if len(extver) == 1:
                #     extname = "COMB"

                # update the master key to be associated with the correct
                # extension
                # spec2D_key[(header["EXTVER"], extname)]["ext"] = ext
                spec2D_key[(extname)]["ext"] = ext

                ## width and height:
                spec2D_key[(extname)]["x_pix"] = naxis1
                spec2D_key[(extname)]["y_pix"] = naxis2

    return spec2D_key


def display_image_in_DS9(fits_filename, frame=None, zoom='to fit', scale='zscale', 
                        cmap='grey', lock_scale=False, regions=None):
    """
    Display a FITS image in DS9 using SAMP instead of XPA.
    
    Parameters:
    -----------
    fits_filename : str or Path
        Path to FITS file to display
    frame : int, optional
        DS9 frame number (SAMP doesn't directly control frames)
    zoom : str or float, optional
        Zoom level ('to fit', 'to width', 'to height', or numeric value)
    scale : str, optional
        Image scaling ('linear', 'log', 'sqrt', 'squared', 'asinh', 'sinh', 'histequ', 'zscale', 'zmax')
    cmap : str, optional
        Colormap name
    lock_scale : bool, optional
        Whether to lock scale across frames
    regions : str, optional
        Path to region file to load
    
    Returns:
    --------
    bool : True if successful, False otherwise
    
    Examples:
    ---------
    # Basic usage
    display_image_in_DS9('image.fits')
    
    # With custom scaling and zoom
    display_image_in_DS9('image.fits', zoom=2, scale='log', cmap='heat')
    
    # With regions
    display_image_in_DS9('image.fits', regions='sources.reg')
    """
    ds9_title = "POPPIES_DIRECT"
    
    try:
        # Get SAMP helper
        samp = get_ds9_samp_helper()

        # Check if file exists
        fits_path = Path(fits_filename)

        if not fits_path.exists():
            logger.error(f"FITS file not found: {fits_filename}")
            return False
        
        # Load the image
        image_name = f"{ds9_title}: {fits_path.name}"
        success = samp.send_image_to_ds9(fits_path, image_name, ds9_target = ds9_title)
        
        if not success:
            return False
        
        # Apply display settings
        # Note: Small delay to ensure image is loaded before applying settings
        time.sleep(0.5)
        
        # Set zoom
        if zoom:
            if isinstance(zoom, (int, float)):
                samp.send_command_to_ds9(f"zoom {zoom}")
            else:
                samp.send_command_to_ds9(f"zoom {zoom}")
        
        # Set scale
        if scale:
            samp.send_command_to_ds9(f"scale {scale}")
        
        # Set colormap
        if cmap:
            samp.send_command_to_ds9(f"cmap {cmap}")
        
        # Lock scale if requested
        if lock_scale:
            samp.send_command_to_ds9("scale lock yes")
        
        # Load regions if provided
        if regions and Path(regions).exists():
            # Note: Region loading via SAMP is not standard, may not work
            region_path = Path(regions).resolve()
            region_params = {
                "url": region_path.as_uri(),
                "name": f"{ds9_title}: {region_path.name}"
            }
            region_message = {
                "samp.mtype": "ds9.regions.load",
                # "samp.mtype": "regions",
                "samp.params": region_params
            }
            try:
                samp.client.notify_all(region_message)
                logger.info(f"Sent region file: {regions}")
            except Exception as e:
                logger.warning(f"Could not load regions via SAMP: {e}")
                logger.info("Consider loading regions manually in DS9")
        
        logger.info(f"Successfully displayed {fits_filename} in DS9")
        return True
        
    except Exception as e:
        logger.error(f"Error displaying image in DS9: {e}")
        return False


def display_images_in_DS9(fits_filenames, tile_mode=True, match_scale=True, 
                         zoom='to fit', scale='zscale', cmap='grey', 
                         delay=0.5, regions=None):
    """
    Display multiple FITS images in DS9 using SAMP.
    
    Parameters:
    -----------
    fits_filenames : list
        List of FITS file paths
    tile_mode : bool, optional
        Whether to use tile mode (arrange images in grid)
    match_scale : bool, optional
        Whether to match scale/zoom across all images
    zoom : str or float, optional
        Zoom level for all images
    scale : str, optional
        Image scaling for all images
    cmap : str, optional
        Colormap for all images
    delay : float, optional
        Delay between loading images (seconds)
    regions : str or list, optional
        Region file(s) to load
    
    Returns:
    --------
    bool : True if all images loaded successfully, False otherwise
    
    Examples:
    ---------
    # Display multiple images in tile mode
    files = ['img1.fits', 'img2.fits', 'img3.fits']
    display_images_in_DS9(files, tile_mode=True, match_scale=True)
    
    # Display with custom settings
    display_images_in_DS9(files, zoom=1.5, scale='log', cmap='heat')
    """
    ds9_title = "POPPIES_DIRECT"
    

    try:   
        # Get SAMP helper
        samp = get_ds9_samp_helper()

        try:
                
            #### Identify the correct client ID:
            labels = samp.identify_client_labels()
            for i, (ids, labs) in enumerate(labels.items()):
                for j,name in enumerate(labs):
                    if name == ds9_title:
                        cl = j

            clientID = labels["client_id"][cl]

        except Exception as e:
            print("Client ID not found due to error: ",e)
            return False

        ####  ####

        # # Set title for first DS9 instance and remember it
        # samp.set_ds9_title_and_remember(ds9_title, ds9_index=0)
        
        # Validate input
        if not fits_filenames:
            logger.error("No FITS files provided")
            return False
        
        # Convert to list if single filename
        if isinstance(fits_filenames, (str, Path)):
            fits_filenames = [fits_filenames]
        
        # Set up DS9 for multiple images
        if tile_mode:
            samp.send_command_to_ds9("tile yes",client_id="c1")
            samp.send_command_to_ds9("tile grid",client_id="c1")
        
        if match_scale:
            samp.send_command_to_ds9("scale lock yes",client_id="c1")
            samp.send_command_to_ds9("zoom lock yes",client_id="c1")
        
        success_count = 0
        
        # Load each image
        for i, (filter_name, filter_images) in enumerate(fits_filenames.items(), start=1):
            for j, image_file in enumerate(filter_images, start=1):

            # for i, fits_file in enumerate(fits_filenames):
                try:
                    fits_path = Path(image_file)
                    if not fits_path.exists():
                        logger.warning(f"FITS file not found: {image_file}")
                        continue
                    
                    # Create new frame for each image (DS9-specific command)
                    if i > 0:  # Don't create frame for first image
                        samp.send_command_to_ds9("frame new",client_id=clientID)
                    
                    # Load the image
                    image_name = f"{ds9_title} {i+1}: {fits_path.name}"
                    if samp.send_image_to_ds9(fits_path, image_name,client_id=clientID):
                        success_count += 1
                        
                        # Apply settings to this image
                        time.sleep(delay)  # Allow time for image to load
                        
                        if zoom:
                            if isinstance(zoom, (int, float)):
                                samp.send_command_to_ds9(f"zoom {zoom}",client_id=clientID)
                            else:
                                samp.send_command_to_ds9(f"zoom {zoom}",client_id=clientID)
                        
                        if scale:
                            samp.send_command_to_ds9(f"scale {scale}",client_id=clientID)
                        
                        if cmap:
                            samp.send_command_to_ds9(f"cmap {cmap}",client_id=clientID)
                        
                        logger.info(f"Loaded image {i}/{len(fits_filenames)}: {image_file}")
                    else:
                        logger.error(f"Failed to load image: {image_file}")
                    
                except Exception as e:
                    logger.error(f"Error loading {image_file}: {e}")
                    continue
        
            # Load regions if provided
            if regions:
                # for i, (filter_name, filter_images) in enumerate(fits_filenames.items(), start=1):
                #     for j, image_file in enumerate(filter_images, start=1):

                region_file = regions[filter_name][j - 1]

                    # for region_file in regions:

                # for region_file in region_files:
                if Path(region_file).exists():

                    samp.send_command_to_ds9(f"region {region_file}",client_id=clientID)

                    # Try to load regions (may not work with all DS9 versions)
                    region_path = Path(region_file).resolve()
                    region_params = {
                        "url": region_path.as_uri(),
                        "name": f"{ds9_title} {i+1}: {region_path.name}"
                    }
                    region_message = {
                        "samp.mtype": "ds9.region.load",
                        "samp.params": region_params
                    }
                    # region_message = {
                    #     "samp.mtype": "regions",
                    #     "samp.params": region_params
                    # }

                    ### commented out for now:
                    # try:
                    #     samp.client.notify_all(region_message)
                    #     logger.info(f"Sent region file: {region_file}")
                    # except Exception as e:
                    #     logger.warning(f"Could not load regions via SAMP: {e}")

        # Final tile arrangement if requested
        if tile_mode and success_count > 1:
            time.sleep(0.5)  # Allow images to settle
            samp.send_command_to_ds9("tile",client_id=clientID)
        
        # if pan_mode:
        #     # print("panning to ra={} and dec={}".format(str(ra[0],str(dec[0]))))
        #     panDirect_POPPIES(samp, ra, dec)

        logger.info(f"Successfully loaded {success_count}/{len(fits_filenames)} images")
        return success_count == len(fits_filenames)

    except Exception as e:
        logger.error(f"Error displaying images in DS9: {e}")
        return False
    
