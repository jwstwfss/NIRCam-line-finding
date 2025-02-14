#### This python script is mostly for DS9 GUIs (as a sort of secondary script)

# from collections import OrderedDict
import os
import numpy as np
from astropy.io import fits

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


def display_spec2D_ds9():
    pass


def display_image_in_DS9(frame_number, image_file, region_file, verbose=True):
    """Display an image in ds9 via the xpa system. Should display an image and
    region file if given.:

    Parameters
    ----------
    frame_number: int
        frame number to display the image
    image_file: str
        image filename or path to display.
    region_file: str
        region file or pathj to display
    """
    # image_file path
    # frame number to add image to
    # add ,region file if ggiven
    ds9_title = "POPPIES_DIRECT"
    
    if image_file is None:
        os.system(f"xpaset -p {ds9_title} frame {frame_number}")
        return

    if verbose:
        msg = f"Opening frame {frame_number}: {os.path.basename(image_file)} | Region File: {region_file}"
        print(msg)

    os.system(f"xpaset -p {ds9_title} frame {frame_number}")  # Specify the frame
    os.system(f"xpaset -p {ds9_title} file {image_file}")

    if region_file:
        os.system(f"xpaset -p {ds9_title} regions load {region_file}")
    # os.system(f"xpaset -p {ds9_title} cmap bb")
    return


def display_images_in_DS9(images, region_files, root=None, verbose=True):
    """Main interface to display multi-images (direct and grism data)
    in DS9. This will produce a tiled display of the direct images (for now).

    TODO: Some of the logic I think can be improved by for not it works.

    Parameters
    ----------
    images: dict
        dictionary contains the image files in different filters to display.
    region_files: dict
        same as the images argument but for regions files to display. Default None
    root: str
        the root or base path the set of files. This to be removed
    verbose:
        messages to output to use if needed. Default is True.

    Returns
    -------
    None
    """

    #F for i, (filter_name, filter_images) in enumerate(images.items(), start=1):
    #     for j, image_file in enumerate(filter_images, start=1):
    #         # print(f"(i, j) | ({i}, {j})")
    #         row = (j - 1) // 3 + 1  # Calculate row number
    #         col = (j - 1) % 3 + 1  # Calculate column number
    #         frame_num = (i - 1) * 3 + (row - 1) * 3 + col
    frame_num = 0
    image_dict = {}
    for i, (filter_name, filter_images) in enumerate(images.items(), start=1):
        for j, image_file in enumerate(filter_images, start=1):
            # print(f"(i, j) | ({i}, {j})")
            # row = (j - 1) // 3 + 1  # Calculate row number
            # col = (j - 1) % 3 + 1  # Calculate column number
            # frame_num = (i - 1) * 3 + (row - 1) * 3 + col
            frame_num += 1
            image_dict[frame_num] = {}

            # Build full paths to image and region files
            # this step should be removed doing to much
            # complete paths should be given already.
            if root:
                image_path = os.path.join(root, image_file)
            else:
                image_path = image_file

            region_file = region_files[filter_name][j - 1] if region_files else None
            # print(f"---DEBUG: WHAT IS THE REGION FILE? - {region_file}")
            # print(f"---DEBUG: WHAT IS THE IMAGE FILE? - {image_file}")

            # add logic to handle missing frames to  keep things symmetric
            # when displaying
            # if image_path is None:
            #     os.system(f"xpaset -p ds9 frame {frame_num}")
            #     continue

            #F # Check if the image file exists
            # if not os.path.isfile(image_path):
            #     if verbose:
            #         print(f"Error: Image file {image_path} does not exist.")
            #     continue

                        # Check if the image file exists
            if not os.path.isfile(image_path):
                if verbose:
                    print(f"Error: Image file {image_path} does not exist.")
                image_dict[frame_num]["img"] = None
                image_dict[frame_num]["reg"] = None
                continue

            # Check if the region file exists
            # print(f"---DEBUG: checking region file - {region_file}")
            if region_file is None:
                is_region_file = False
            else:
                is_region_file = os.path.isfile(region_file)

            # print(f"---DEBUG: IS REGION FILE - {is_region_file}")

            if region_file and not is_region_file:
                if verbose:
                    print(f"Error: Region file {region_file} does not exist.")
                region_file = None

            # Display the image
            image_dict[frame_num]["img"] = image_path
            image_dict[frame_num]["reg"] = region_file

    for frame_num in image_dict:
        display_image_in_DS9(frame_num,
                             image_dict[frame_num]["img"],
                             image_dict[frame_num]["reg"], verbose=verbose)
        
            #F # Display the image
            # display_image_in_DS9(frame_num, image_path, region_file, verbose=verbose)

    ds9_title = "POPPIES_DIRECT"
    # Go to frame 1
    os.system(f"xpaset -p {ds9_title} frame 1")

    #F # Configure additional settings
    # os.system(f"xpaset -p {ds9_title} tile")
    # # os.system(f"xpaset -p {ds9_title} lock frame wcs")
    # os.system(f"xpaset -p {ds9_title} lock scale")
    # # os.system(f"xpaset -p {ds9_title} zoom to fit")
    # os.system(f"xpaset -p {ds9_title} scale mode zscale")
    # os.system(f"xpaset -p {ds9_title} lock colorbar")
    # #os.system(f"xpaset -p {ds9_title} lock scalelimits")

    # Configure additional settings
    os.system(f"xpaset -p {ds9_title} tile")
    # os.system(f"xpaset -p {ds9_title} match frame wcs")
    # os.system(f"xpaset -p {ds9_title} lock scale")
    # os.system(f"xpaset -p {ds9_title} zoom to fit")
    os.system(f"xpaset -p {ds9_title} scale mode zscale")
    os.system(f"xpaset -p {ds9_title} lock colorbar")
    #os.system(f"xpaset -p {ds9_title} lock scalelimits")

    # Figured out the angle is ~250 for PAR 28. Need to check if this is true for all.
    ## FH updated 1/24/25:
    # for framename in range(1,10):
    
    for framename in range(1,2):
    
        os.system(f"xpaset -p {ds9_title} frame " + str(framename))
        os.system(f"xpaset -p {ds9_title} wcs sky ecliptic")
        os.system(f"xpaset -p {ds9_title} zoom to fit")


    # os.system(f"xpaset -p {ds9_title} frame 3")
    os.system(f"xpaset -p {ds9_title} match frame image")

    # if I match at open -- does it stay matched going forward?
    # Go to frame 1
    # os.system(f"xpaset -p {ds9_title} frame 3")
    # os.system(f"xpaset -p {ds9_title} match frame wcs")

    return
