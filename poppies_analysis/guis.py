#### This python script is mostly for DS9 GUIs

import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
import os
from poppies_analysis import *
from distutils.sysconfig import *
from array import array
from glob import glob
import pandas as pd
from ast import literal_eval


# POPPIES gui helpers
from guis_helpers import extract_image_extensions_key
from guis_helpers import display_image_in_DS9, display_images_in_DS9

# xpa resources:
# https://www.astro.louisville.edu/software/xmccd/archive/xmccd-4.1/xmccd-4.1e/docs/xpa/xpa.pdf
# http://ds9.si.edu/doc/ref/xpa.html#lock    

### FH updating 1/13/25
def showSpec2D_POPPIES(parno, obid, filter = "F444W", path_to_data = "", zsc = "squared"):
    """Display spec2D cutouts in DS9.

    This can handle missed data in a given filter since there is a key
    to check against given the nature of how the spec2D files are setup.
    Frames will be made for missing data but will show nothing to indicate
    that the data is missing and that display layout is preserved.
    """

    # try:
    #     secats = glob(path_to_data + "/Par" + str(parno) + "/DATA/DIRECT_GRISM/Par*phot*.fits")  # MDR 2022/05/17
    # except:
    #     secats = glob(path_to_data + "/Par" + str(parno) + "/Products/Par*phot*.fits")  # KVN allowing for different path structure (?)


    def parse_filename(path_to_data, parno, obid):
        # return path_to_data + f"/Par{parno}/spec2D/Par{parno}_{obid:05d}.2D.fits"
        f1 = glob(path_to_data + f"/Par{parno}/Spec1D2D/A_{filter}_R_{obid}.V4.fits")  
        f2 = glob(path_to_data + f"/Par{parno}/Spec1D2D/A_{filter}_C_{obid}.V4.fits") 

        f3 = glob(path_to_data + f"/Par{parno}/Spec1D2D/B_{filter}_R_{obid}.V4.fits")  
        f4 = glob(path_to_data + f"/Par{parno}/Spec1D2D/B_{filter}_C_{obid}.V4.fits")         
        
        return f1,f2,f3,f4


    spec2D_files = parse_filename(path_to_data, parno, obid)

    #check if given grism 1d files exist or not:
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

    ### first for GRISM-R:
    # spec2D_file = parse_filename(path_to_data, parno, obid)[0][0]

    SPEC2D_TITLE_DS9 = "POPPIES_spec2D"


    if grismr_file is not None:
            
        # for a given data file extract the data information that we want to be displayed in ds9
        spec2D_key_DS9 = extract_image_extensions_key(grismr_file)


        # display in ds9 instance of a given "title"

        # now display the results.
        for key, _ in spec2D_key_DS9.items():

            # get key-values for frame id and extension
            frame_id = int(spec2D_key_DS9[key]["frame_id"]) - 1    ## -1 bc first extension is 1D spec
            ext = spec2D_key_DS9[key]["ext"]

            # # display frame
            if ext is not None:
                
                #width and height of 2D spec trace
                xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]

                xpix,ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2*ypix_raw))

                # frame id
                command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}"

                os.system(command)

                command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

                os.system(command)

                command = f"xpaset -p {SPEC2D_TITLE_DS9} fits {grismr_file}[{ext}]"

                os.system(command)

                os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
                os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")
                
        # else:
        #     command = f"xpaset -p {SPEC2D_TITLE_DS9} frame clear"
        #     os.system(command)
    
    else:
        for frame_id in [1,2,3]:
            # frame id
            os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}")
            os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame clear")

    # # formating display into grid
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile")
    # # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid mode manual")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile mode row")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid layout 1 6")

    # # config display properties
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 1")
    # os.system("xpaset -p {SPEC2D_TITLE_DS9} lock frame physical")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock colorbar")

    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode user")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale yes")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits yes")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap grey")


    # for fno in [1,2,3]:
    #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " +str(fno))
    #     # OLD: 39.5 and 32.5 are just 1/2 cutout length & 1/2 cutout width. Will need to be changed for NIRCam
    #     # The length is set to the length of the cutout -8 (for 4 pix on each end)
    #     # For NIRCam F444W, image size is 1408 x 200
    #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 710 105 100 10# color=green} ")        
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale squared")
    #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale")
    #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits yes")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
        # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1.5")

        # command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

        # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh")

    ### Now for GRISM-C:
    # spec2D_file = parse_filename(path_to_data, parno, obid)[1][0]


    if grismc_file is not None:

        # for a given data file extract the data information that we want to be displayed in ds9
        spec2D_key_DS9 = extract_image_extensions_key(grismc_file)

        # now display the results.
        for key, _ in spec2D_key_DS9.items():

            # get key-values for frame id and extension
            frame_id = (int(spec2D_key_DS9[key]["frame_id"]) - 1) + 3    ## +3 bc new oreintation here
            ext = spec2D_key_DS9[key]["ext"]

            # # display frame
            if ext is not None:
                
                # print('GC ', frame_id)

                #width and height of 2D spec trace
                xpix_raw, ypix_raw = spec2D_key_DS9[key]["x_pix"], spec2D_key_DS9[key]["y_pix"]

                xpix,ypix = int(xpix_raw + (0.2 * xpix_raw)), int(ypix_raw + (2*ypix_raw))

                # frame id
                command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}"

                os.system(command)

                command = f"xpaset -p {SPEC2D_TITLE_DS9} frame refresh"

                os.system(command)

                command = f"xpaset -p {SPEC2D_TITLE_DS9} fits {grismc_file}[{ext}]"

                os.system(command)

                os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")
                os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")



            # else:
            #     command = f"xpaset -p {SPEC2D_TITLE_DS9} frame clear"
            #     os.system(command)


    else:
        for frame_id in [4,5,6]:
            # frame id
            os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}")
            os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame clear")

    # formating display into grid
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid mode manual")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile mode row")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid layout 1 6")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width 1600")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height 1200")

    try:  #FH modified 3/6/25
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")
    except Exception as e:
        print('Could not correctly display spec2D, Reason: ', e)
        pass

    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale {zsc}")

    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock colorbar")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")

    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")

    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 7")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")

    # # for fno in np.arange(6,0,-1):

    for fno in np.arange(1,8,1):

        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " + str(fno))
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame center")

        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom 1")
        
        # FH 2/17/25: display wavelength beside alpha
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} wcs fk5")

        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} wcs skyformat degrees")


    #     # For NIRCam F444W, image size is 1408 x 200
    #     # os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 710 105 100 10# color=green} ")        
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} width {xpix}")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} height {ypix}")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode zscale")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale mode 99.5")
    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} scale {zsc}")

    #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap heat")

        # ## Need a 180 degree rotation for R-spectra:
        # if fno in [1,2,3]:
        #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} rotate 180")

        # if fno==7:
        #     os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")



    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 7")
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame delete")

    # # Go to the first frame:    
    # os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 1")


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

# FH updated 2/7/25
def showDirect_POPPIES(parno, filter="F444W", path_to_data=""):
    """Displays direct images for each of the filters"""

    ### KVN's quick fix to images having different names in different fields
    # specify the images to be displayed in DS9
    # grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*gr150*_drz_sci.fits')
    ### FH updated 1/6/25
    # grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*A_F444W_i2d.fits')
    # grism_file = str(grism_file[0]).split(path_to_data)[1]
    # grism_file_2 = str(grism_file).split('DATA/')[1]


    image_files = glob(path_to_data + 'Par'+str(parno)+'/DirectImages/*_{}_i2d.fits'.format(filter))
    
    image_file = image_files[0]

    image_file_ext = str(image_file).split('DirectImages/')
    image_file_2 = image_file_ext[1]

    images = {
        f"{filter}": [
            f"{image_file_2}"
            # f"Par{parno}_{grism_file_ext}_f115w-gr150c_drz_sci.fits",
            # f"Par{parno}_{grism_file_ext}_f115w-gr150r_drz_sci.fits",
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
                print(f"Warning: File {filename} not found; assuming it lives in path_to_data/data/Par#/DirectImages.")
                full_path = find_file(path_to_data + '/Par'+str(parno)+'/DirectImages/', filename)
                paths.append(full_path)
                print(paths)
                
        image_paths[filter_name] = paths

    
    # specify the region files to show along with images, if any.
    # region_files = {
    #     "f115w": [f"Par{parno}regions_phot.reg", f"Par{parno}F115c_grism.reg", f"Par{parno}F115r_grism.reg"],
    #     "f150w": [f"Par{parno}regions_phot.reg", f"Par{parno}F150c_grism.reg", f"Par{parno}F150r_grism.reg"],
    #     "f200w": [f"Par{parno}regions_phot.reg", f"Par{parno}F200c_grism.reg", f"Par{parno}F200r_grism.reg"],
    # }
    region_files = {
        f"{filter}": [f"Par{parno}regions_phot_{filter}.reg"]
    }

    # find the full paths of the region files
    region_paths = {}
    for filter_name, filenames in region_files.items():
        paths = []
        for filename in filenames:
            full_path = find_file(path_to_data, filename)
            # print(full_path)
            if full_path:
                paths.append(full_path)
            else:
                # KVN: assuming same pathing structure
                print(f"Warning: File {filename} not found; assuming it lives in path_to_data/data/Par#/DirectImages.")
                full_path = find_file(path_to_data + '/Par'+str(parno)+'/DirectImages/', filename)
                paths.append(full_path)
                
        region_paths[filter_name] = paths

    display_images_in_DS9(image_paths, region_files=region_paths)

    return


## FH updated 1/24/25
def panDirect_POPPIES(ra, dec):
    # print(ra, dec)

    ds9_title = "POPPIES_DIRECT"
    # for fno in [1,4,7]:

    ## will do loop when more filters are added
    # for fno in [1]:
    cmd = f"xpaset -p {ds9_title} frame 1"  # + str(fno)
    os.system(cmd)
    cmd = f"xpaset -p {ds9_title} pan to {ra[0]} {dec[0]} wcs degrees"
    os.system(cmd)

    # # zoom
    # cmd = f"xpaset -p {ds9_title} zoom 2"
    # os.system(cmd)

    # os.system(f"xpaset -p {ds9_title} zoom in")
    os.system(f"xpaset -p {ds9_title} zoom to 4")


# written by KVN June 2024
# First pull info from region file
def getRegionFileInfo(parno, filt, path_to_data):

        if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+str(filt)+'_grism.reg'):
                df = pd.read_csv(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+str(filt)+'_grism.reg', sep=' ', engine='python', header=None)
        
                # pull the id and x & y coord from the region file
                obj_ids = [i for i in df[3].apply(lambda z: int(z.split('=')[-1].strip('}{')))]
                obj_xs = df[0].apply(lambda z: literal_eval(z.lstrip('box'))[0])
                obj_ys = df[0].apply(lambda z: literal_eval(z.lstrip('box'))[1])
        else:
                print('No region file found for panning to object')
                obj_ids = None; obj_xs = None; obj_ys = None

        return obj_ids, obj_xs, obj_ys


# written by KVN June 2024
# Second, use info from region file to pan to each grism
# FH updated 1/7/25
def panDispersed_POPPIES(objid, parno, path_to_data):

    # grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*gr150*_drz_sci.fits')
    # grism_file = str(grism_file[0]).split(path_to_data)[1]
    # grism_file_ext = str(grism_file).split('_')[1]

    grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*A_F444W_i2d.fits')
    grism_file = str(grism_file[0]).split(path_to_data)[1]
    # grism_file_ext = str(grism_file).split('_')[1]
    grism_file_2 = str(grism_file).split('DATA/')[1]

    ds9_title = "POPPIES_DIRECT"
    # since tiles are already assigned, setup by tile rather than grism
    if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/'+str(grism_file_2)):
        fno = 2
        obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F444c')
        
        # specify the frame
        cmd = f"xpaset -p {ds9_title} frame " + str(fno)
        os.system(cmd)

        # pan to the coords from the region file for obj with objid
        get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
        cmd = "xpaset -p POPPIES_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
        os.system(cmd)
        cmd = "xpaset -p POPPIES_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
        os.system(cmd)


    if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/'+str(grism_file_2)):
        fno = 3
        obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F444r')

        # specify the frame
        cmd = f"xpaset -p {ds9_title} frame " + str(fno)
        os.system(cmd)

        # pan to the coords from the region file for obj with objid
        get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
        cmd = "xpaset -p POPPIES_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
        os.system(cmd)


    # if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150c_drz_sci.fits'):
    #     fno = 5
    #     obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F150c')

    #     # specify the frame
    #     cmd = f"xpaset -p {ds9_title} frame " + str(fno)
    #     os.system(cmd)

    #     # pan to the coords from the region file for obj with objid
    #     get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
    #     cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
    #     os.system(cmd)
    
    # if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150r_drz_sci.fits'):
    #     fno = 6
    #     obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F150r')

    #     # specify the frame
    #     cmd = f"xpaset -p {ds9_title} frame " + str(fno)
    #     os.system(cmd)

    #     # pan to the coords from the region file for obj with objid
    #     get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
    #     cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
    #     os.system(cmd)

    # if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150c_drz_sci.fits'):
    #     fno = 8
    #     obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F200c')

    #     # specify the frame
    #     cmd = f"xpaset -p {ds9_title} frame " + str(fno)
    #     os.system(cmd)

    #     # pan to the coords from the region file for obj with objid
    #     get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
    #     cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
    #     os.system(cmd)

    # if os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150r_drz_sci.fits'):
    #     fno = 9
    #     obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_data=path_to_data, filt='F200r')

    #     # specify the frame
    #     cmd = f"xpaset -p {ds9_title} frame " + str(fno)
    #     os.system(cmd)

    #     # pan to the coords from the region file for obj with objid
    #     get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
    #     cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
    #     os.system(cmd)



# written by KVN June 2024
# DO NOT use this one. Most of below was testing. Faster/Better to get the coords from the region file
# Not the best way to do this... but the positions are in RA/DEC
# so I convert from RA/DEC to pixels then back to RA/DEC
def panDispersed_PASSAGE_doNOTuse(ra, dec, parno, path_to_wisp=""):

    grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*gr150*_drz_sci.fits')
    grism_file_ext = grism_file[0].split('_')[1]

    pix_per_um = 1 / (1e-4 * 46.934)  # GR150R 47.015 for GR150C
    tweak = 7

    ds9_title = "PASSAGE_DIRECT"    
    # since tiles are already assigned, setup by tile rather than grism
    for fno in [2,3]:
        if fno == 2 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f115: ', xx,yy)

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0# -6 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        elif fno == 3 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150r_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150r_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position

            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset,  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)

    for fno in [5,6]:
        if fno == 5 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f150: ', xx,yy)

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0 #55 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        else: # R orientation
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset,  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)


    for fno in [8,9]:
        if fno == 8 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f200: ',xx,yy)
            
            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0#147 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)

        elif fno == 9 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150r_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150r_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0#147 + 0.135 * pix_per_um + tweak
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        
### FH updated 12/24/24
def show2dNEW(
    grism,
    parno,
    obid,
    zeroarr,
    user,
    trans,
    zran1=None,
    zran2=None,
    path_to_data=" ",
):
    # In version 1.0, will first look for wavelength-calibrated stamps in the G1??_DRIZZLE directories; failing this, will default to old stamps
    # zero and first order positions
    #    firstx = firstarr['x']
    #    firsty = firstarr['y']
    #    firstlen = firstarr['len']
    #    firstwid = firstarr['width']
    #    firstid = firstarr['objid']
    zerox = zeroarr["x"]
    zeroy = zeroarr["y"]
    zeroid = zeroarr["objid"]
    zerora = zeroarr["ra"]
    zerodec = zeroarr["dec"]
    zeromag = zeroarr["mag_auto"]

    dims = ()
    zrad = 10.0
    workingdir = os.getcwd()
    par_root_dir = "/"
    if path_to_data == " ":
        dirpts = workingdir.split("/")[1:-1]
        for pdir in dirpts:
            par_root_dir = par_root_dir + pdir + "/"
        # path2dl=par_root_dir + grism + '_DRIZZLE/aXeWFC3_' +grism + '_mef_ID'+str(obid)+'.fits'
        path2dl = (
            par_root_dir
            + "spec1D/"
            + "A_{}_GRISMR_{}".format(str(grism),str(obid))
            + ".V4.fits"
        )
    else:
        # path2dl = path_to_data + '/Par' + str(parno) +'/' + grism + '_DRIZZLE/aXeWFC3_' +grism + '_mef_ID'+str(obid)+'.fits'
        path2dl = (
            path_to_data
            + "Par{}".format(str(parno))
            + "/spec1D/"
            + "A_{}_GRISMR_{}".format(str(grism),str(obid))
            + ".V4.fits"
        )
    
    if os.path.exists(path2dl) == 1:
        path2d = path2dl
    else:

        if path_to_data == " ":
            path2d = (
            par_root_dir
            + "/spec1D/"
            + "A_{}_GRISMR_{}".format(str(grism),str(obid))
            + ".V4.fits"
            )
        else:

            # path2d = path_to_data + '/Par' + str(parno) + '/Stamps/Par' + str(parno) + '_'+grism+'_BEAM_'+str(obid)+'A.fits'
            path2d = (
            path_to_data
            + "/Par{}".format(str(parno))
            + "/spec1D/"
            + "A_{}_GRISMR_{}".format(str(grism),str(obid))
            + ".V4.fits"
            )

    pix_per_um = 1 / (4e-4 * 10)  # Needs changing for NIRCam: GR150R 47.015 for GR150C

    tweak = 7
    tweak_zp = 7
    if grism == "F200W":
        frameno = "3"
        maglimit = 26.0
        input_grism = (
            path_to_data
            + "Par"
            + str(parno)
            + "/DATA/Par"
            + str(parno)
            + "_f200w-gr150r_drz_sci.fits"
        )
        zp_offset = -210
        obj_offset = (
            147 + 0.238 * pix_per_um + tweak
        )  # First number is pixels to blue edge (poorly defined). Second num is wave distance to center in um.
    elif grism == "F150W":
        frameno = "2"
        maglimit = 24.5
        input_grism = (
            path_to_data
            + "Par"
            + str(parno)
            + "/DATA/Par"
            + str(parno)
            + "_f150w-gr150r_drz_sci.fits"
        )
        zp_offset = -219 + tweak_zp
        obj_offset = 55 + 0.171 * pix_per_um + tweak
    elif grism == "F444W":
        frameno = "1"
        maglimit = 23.5
        input_grism = (
        path_to_data + 
        '/Par'+str(parno)+
        '/DATA/*A_F444W_i2d.fits'
        )
        zp_offset = -225.5 + tweak_zp
        obj_offset = -6 + 0.135 * pix_per_um + tweak

    direct1 = (
        path_to_data + 
        '/Par'+str(parno)+
        '/DATA/*A_F444W_i2d.fits'
    )
    direct2 = (
        path_to_data
        + "Par"
        + str(parno)
        + "/DATA/"
        + "Par"
        + str(parno)
        + "_f150w_drz_sci.fits"
    )
    direct3 = (
        path_to_data
        + "Par"
        + str(parno)
        + "/DATA/"
        + "Par"
        + str(parno)
        + "_f200w_drz_sci.fits"
    )

    if os.path.exists(direct1) == 1:
        input_direct = direct1
    elif os.path.exists(direct2) == 1:
        input_direct = direct2
    elif os.path.exists(direct3) == 1:
        input_direct = direct3

    # tiltedgrism = 0
    # extension_check = 0

    ## FH commenting these out for now
    # if os.path.exists(path2d) == 1:
    #     infits = fits.open(path2d)
    #     dal = len(infits)
    #     doug = list(range(1, dal))
    #     for tracker in doug:
    #         if grism == infits[tracker].header["EXTVER"]:
    #             extension_check = 1

    # if extension_check:
    #     cutout_index = infits.index_of(("SCI", grism))
    #     grinfits = fits.open(input_grism)
    #     hdr = grinfits[0].header
    #     if "PA_APER" in hdr:
    #         pa_aper = hdr["PA_APER"]
    #     else:
    #         newfits = fits.open(input_direct)
    #         nhdr = newfits[0].header
    #         pa_aper = nhdr["PA_APER"]
    #         tiltedgrism = 1
    #     # print("PA_APER: ", pa_aper)
    #     wcs_in = WCS(header=grinfits[0].header)
    #     ra = array("d", zerora)
    #     dec = array("d", zerodec)
    #     x, y = wcs_in.all_world2pix(ra, dec, 0)
    #     x_obj = x[obid - 1]
    #     y_obj = y[obid - 1]

    #     # This section probably needs a grism orientation check G150R vs G150C
    #     if tiltedgrism:
    #         print(
    #             "Attention: Reference grism has been tilted from original orientation."
    #         )
    #         pa_rad = -pa_aper / (180 / math.pi)
    #         #  print(pa_rad,zp_offset*math.cos(pa_rad),zp_offset*math.sin(pa_rad))
    #         y_zp = y - zp_offset * math.cos(pa_rad)
    #         x_zp = x - zp_offset * math.sin(pa_rad)
    #         x_spec = x_obj - obj_offset * math.sin(pa_rad)
    #         y_spec = y_obj - obj_offset * math.cos(pa_rad)
    #         x_test = x - obj_offset * math.sin(pa_rad)
    #         y_test = y - obj_offset * math.cos(pa_rad)
    #     else:
    #         y_zp = y - zp_offset
    #         x_zp = x
    #         x_spec = x_obj
    #         y_spec = y_obj - obj_offset
    #         x_test = x
    #         y_test = y - obj_offset

    #     ### changing to read in 1st data extension ###
    #     # darr=infits[-1].data
    #     hdr = infits[cutout_index].header
    #     darr = infits[cutout_index].data
    #     rms = np.std(darr)
    #     dims = darr.shape

    #     y_c = 0.5 * (dims[0] + 1)
    #     x_c = 0.5 * (dims[1] + 1)
    #     x_delta = y_c - x_spec
    #     y_delta = x_c - y_spec
    #     x_zp_new = dims[1] - y_zp - y_delta
    #     y_zp_new = x_delta + x_zp
    #     x_test_new = dims[1] - y_test - y_delta
    #     y_test_new = x_delta + x_test
    #     #     print(x_test[obid-1],y_test[obid-1],x_test_new[obid-1],y_test_new[obid-1],dims[1],dims[0])
    #     xmax = dims[1]
    #     ymax = dims[0]
    #     infits.close()
    #     grinfits.close()

        # if zran1 == None:
        #     zran1 = -1 * rms
        # if zran2 == None:
        #     zran2 = 5 * rms
    
    if (os.path.exists(path2d) == 0):

    # elif (os.path.exists(path2d) == 0) or (extension_check == 0):
        print("%s stamp not found." % (grism))
        return False

    ###### Commenting out everything below for now ######
    # COMMENTING OUT ZERO ORDER DISPLAY -START

    ### USING THE DRIZZLE TRANSFORMATIONS TO GET ZEROTH ORDERS ###
    #    _cx = np.array([xcoo for xcoo in zerox]) - hdr['BB0X'] - 1
    #    _cy = np.array([ycoo for ycoo in zeroy]) - hdr['BB0Y'] - 1
    #    cx = hdr['D001OUXC'] + (hdr['DRZ00'] + hdr['DRZ01']*(_cx-hdr['D001INXC']) + hdr['DRZ02']*(_cy-hdr['D001INYC']))
    #    cy = hdr['D001OUYC'] + (hdr['DRZ10'] + hdr['DRZ11']*(_cx-hdr['D001INXC']) + hdr['DRZ12']*(_cy-hdr['D001INYC']))
    # convert to (Angs,arcsec) coords
    #    cx = (cx - hdr['CRPIX1'])*hdr['CDELT1'] + hdr['CRVAL1']
    #    cy = (cy - hdr['CRPIX2'])*hdr['CDELT2'] + hdr['CRVAL2']
    #    rad = 5 * hdr['CD1_1']

    # if path_to_data != " ":
    #     par_root_dir = path_to_data + "/Par" + str(parno) + "/"
    # outcoo = (
    #     par_root_dir + "Spectra/temp_zero_coords_%s_" % user + str(frameno) + ".reg"
    # )

    # if os.path.exists(outcoo) == 1:
    #     os.unlink(outcoo)
    # f = open(outcoo, "w")
    # # f.write('wcs;\n')
    # xmax1 = 3100
    # ymax1 = 3100

    # test = np.where(
    #     (x_zp_new > 0)
    #     & (x_zp_new < xmax)
    #     & (y_zp_new > 0)
    #     & (y_zp_new < ymax)
    #     & (zeromag < maglimit)
    # )
    # #  test=np.where(((x_zp > 0) & (x_zp < xmax1) & (y_zp > 0) & (y_zp < ymax1) & (zeromag < maglimit)))
    # #    test=np.where(((x > 0) & (x < xmax1) & (y > 0) & (y < ymax1) & (zeromag < maglimit)))
    # cross_size = 50
    # rad = 5
    # f.write(
    #     "point(%.2f,%.2f) # point=cross %i color=yellow text={%.1f}\n"
    #     % (x_test_new[obid - 1], y_test_new[obid - 1], cross_size, zeromag[obid - 1])
    # )
    # if test[0] != []:
    #     for j in range(0, len(test[0])):
    #         inputx = x_zp_new[test[0][j]]
    #         #    inputx=x_zp[test[0][j]]
    #         #    inputx=x[test[0][j]]
    #         inputy = y_zp_new[test[0][j]]
    #         #     inputy=y_zp[test[0][j]]
    #         #     inputy=y[test[0][j]]
    #         inputID = zeroid[test[0][j]]
    #         inputmag = zeromag[test[0][j]]
    #         #      testx=x_test[test[0][j]]
    #         #      testy=y_test[test[0][j]]
    #         print(inputx, inputy, rad, inputID)
    #         f.write(
    #             "circle(%.2f,%.2f,%.1f) # color=red text={%.1f;%s}\n"
    #             % (inputx, inputy, rad, inputmag, inputID)
    #         )
    # #      f.write('circle(%.2f,%.2f,%.1f) # color=blue text={%.1f;%s}\n' % (testx,testy,rad,inputmag,inputID))
    # f.close()

    # COMMENTING OUT ZERO ORDER DISPLAY -END

    ### THIS WAS ALL FOR THE OLD TRANSFORMATION ###
    #    matchind=0
    #    i=0
    #    for fid in firstid:
    #        if fid==obid:
    #            matchind=i
    #            break
    #        i=i+1
    #    xmin=firstx[matchind]-firstlen[matchind]/2.0
    #    xmax=firstx[matchind]+firstlen[matchind]/2.0
    #    ymin=firsty[matchind]-firstwid[matchind]/2.0
    #    ymax=firsty[matchind]+firstwid[matchind]/2.0
    #
    #    numzer=0
    #    if len(dims)>0:
    #        xdim=float(max(dims))
    #        ydim=float(min(dims))
    #    outcoo=par_root_dir+"Spectra/temp_zero_coords.reg"
    #    if os.path.exists(outcoo)==1:
    #        os.unlink(outcoo)
    #    coordout=open(outcoo,'w')
    #    print >>coordout, "image"
    #    for j in range(len(zerox)):
    #        # MB: use dims of stamp, rather than size of 1st order in region file
    #        if zerox[j] >= (firstx[matchind] - xdim/2.) and \
    #           zerox[j] <= (firstx[matchind] + xdim/2.) and \
    #           zeroy[j] >= (firsty[matchind] - ydim/2.) and \
    #           zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G102':
    ##    	if zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G102':
    #    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+212./2.0-13,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+3.6,zeroid[j])
    #
    #        elif zerox[j] >= (firstx[matchind] - xdim/2.) and \
    #             zerox[j] <= (firstx[matchind] + xdim/2.) and \
    #             zeroy[j] >= (firsty[matchind] - ydim/2.) and \
    #             zeroy[j] <= (firsty[matchind] + ydim/2.) and grism=='G141':
    ##    	elif  zerox[j]>=(xmin-zrad) and zerox[j]<=(xmax+zrad) and zeroy[j]>=(ymin-zrad) and zeroy[j]<=(ymax+zrad) and grism=='G141':
    #    		print >>coordout, "circle(%.2f,%.2f,5.0) # text={%s}" % (zerox[j]/1.7-firstx[matchind]/1.7+184./2.0,zeroy[j]/1.6-firsty[matchind]/1.6+ydim/2.0+0.6,zeroid[j])
    #    numzer=numzer+1
    #    coordout.close()

    # if trans == "log":
    #     zscale = "log"
    # else:
    #     zscale = "linear"
    # cmd = "xpaset -p ds9 frame " + frameno
    # os.system(cmd)
    # cmd = "xpaset -p ds9 file " + path2d + "[" + str(cutout_index) + "]"
    # os.system(cmd)
    # cmd = "xpaset -p ds9 scale limits " + str(zran1) + " " + str(zran2)
    # os.system(cmd)
    # cmd = "xpaset -p ds9 scale " + zscale
    # os.system(cmd)

    # cmd='xpaset -p ds9 zoom to fit'
    # os.system(cmd)
    # HAVE COMMENTED OUT overplotting of some REGIONS files

    #  if test[0] != []:
    # cmd = (
    #     "xpaset -p ds9 regions file " + outcoo
    # )  # par_root_dir+ 'Spectra/temp_zero_coords_%s.reg'%user
    # os.system(cmd)

    # MR

    #  if frameno=='1':
    #     cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G102_trace.reg'
    #  if frameno=='2':
    #     cmd='xpaset -p ds9 regions file '+par_root_dir+ 'Spectra/G141_trace.reg'
    #  os.system(cmd)

    ###### ######

    return



# PASSAGE version:
def showDirectNEW(obid, parno, g102zeroarr, load_image=False, path_to_data=" "):
    """
    Removed lineno, which was only used to check whether the images
    should be reloaded.

    """
    obid = int(obid)
    workingdir = os.getcwd()
    dirpts = workingdir.split("/")[1:-1]
    par_root_dir = "/"
    if path_to_data == " ":
        for pdir in dirpts:
            par_root_dir = par_root_dir + pdir + "/"
    else:
        par_root_dir = path_to_data + "Par" + str(parno) + "/"

    path2direct = par_root_dir + "DATA/"
    path110 = path2direct + "Par" + str(parno) + "_f115w_drz_sci.fits"
    path140 = path2direct + "Par" + str(parno) + "_f150w_drz_sci.fits"
    path160 = path2direct + "Par" + str(parno) + "_f200w_drz_sci.fits"
    path140cat = path2direct + "DIRECT_GRISM/" + "Par" + str(parno) + "photcat.fits"
    path160cat = path2direct + "fin_F160.cat"
    path110cat = path2direct + "fin_F110.cat"

    # print(obid)
    inputx = g102zeroarr["x"]
    inputy = g102zeroarr["y"]

    if (
        os.path.exists(path110) == 0
        and os.path.exists(path140) == 0
        and os.path.exists(path160) == 0
    ):
        print("No Direct Images Found.")
        return 0
    # For PASSAGE have just the one object catalog. Feeding it in from measure_z_interactive, not read in here.
    # Cutting this section

    #        if os.path.exists(path140cat)==1:
    #        infHcat=open(path140cat,'r')
    #    elif os.path.exists(path160cat)==1:
    #        infHcat=open(path160cat,'r')
    #    elif os.path.exists(path110cat)==1:
    # if no Hband, use F110 for catalog
    #        infHcat=open(path110cat,'r')
    #    else:
    #        return 0
    xcen, ycen = -1, -1
    #   for line in infHcat:
    #       if line[0]!='#':
    #           entries=line.split()
    #          if int(entries[1])==obid:
    #             xcenter,ycenter=float(entries[7]),float(entries[8])
    #             hexcoo=[entries[7],entries[8]]
    # infHcat.close()
    hexcoo = [str(inputx[obid - 1]), str(inputy[obid - 1])]

    # Added getting PA angle of each image.

    '''
    # load the direct images
    if load_image:
        print("IMAGE IS BEING RELOADED AND ROTATION APPLIED!")
        if os.path.exists(path110) == 1:
            hdu110 = fits.open(path110)
            hdr = hdu110[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 4"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path110
            os.system(cmd)
            ### using F110_drz.reg with F110W_drz.fits
            #    cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg'
            #    os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            os.system(cmd)
        if os.path.exists(path140) == 1:
            hdu140 = fits.open(path140)
            hdr = hdu140[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 5"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path140
            os.system(cmd)
            #      cmd='xpaset -p ds9 regions file '+par_root_dir+ 'DATA/DIRECT_GRISM/F140_drz.reg'
            #     os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            os.system(cmd)
        if os.path.exists(path160) == 1:  # was elif
            hdu160 = fits.open(path160)
            hdr = hdu160[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 6"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path160
            print(cmd)
            os.system(cmd)
            #      cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg'
            #     os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            print(cmd)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            print(cmd)
            os.system(cmd)
    '''
    # pan to the coordinates of this object
    if os.path.exists(path110):
        panDirect(hexcoo[0], hexcoo[1], grism="F115W")
    panDirect(hexcoo[0], hexcoo[1])
    if os.path.exists(path140):
        panDirect(hexcoo[0], hexcoo[1], grism="F150W")
    panDirect(hexcoo[0], hexcoo[1])
    if os.path.exists(path160):
        panDirect(hexcoo[0], hexcoo[1], grism="F200W")
    panDirect(hexcoo[0], hexcoo[1])


def showDispersed(obid, parno, load_image=False, path_to_data=" "):  # MB
    """
    Removed lineno, which was only used to check whether the images
    should be reloaded.
    """
    obid = int(obid)
    workingdir = os.getcwd()
    dirpts = workingdir.split("/")[1:-1]
    par_root_dir = "/"
    # for pdir in dirpts:
    #    par_root_dir= par_root_dir +pdir + '/'

    if path_to_data == " ":
        for pdir in dirpts:
            par_root_dir = par_root_dir + pdir + "/"
    else:
        par_root_dir = path_to_data + "/Par" + str(parno) + "/"

    path2dispersed = par_root_dir + "Spectra/DATA/DIRECT_GRISM/"
    ### Using G102.fits instead of G102_drz.fits ###
    path102 = (
        path2dispersed + "Par" + str(parno) + "_" + str(obid).zfill(5) + ".2D.fits"
    )
    path141 = (
        path2dispersed + "Par" + str(parno) + "_" + str(obid).zfill(5) + ".2D.fits"
    )
    path102_0reg = os.path.join(path2dispersed, "G102_0th.reg")
    path102_1reg = os.path.join(path2dispersed, "G102_1st.reg")
    path141_0reg = os.path.join(path2dispersed, "G141_0th.reg")
    path141_1reg = os.path.join(path2dispersed, "G141_1st.reg")
    if os.path.exists(path102) == 0 and os.path.exists(path141) == 0:
        print("Looking for objects here: ", path102)
        print("No Grism Images Found.")
        return 0
    # get center of 1st order
    ### Using same syntax as getzeroorders for consistency ###
    print(path102_1reg)
    if os.path.exists(path102_1reg) == 1:
        reg102 = open(path102_1reg, "r")
        x102, y102 = -1, -1
        for line in reg102:
            # using same syntax as getzeroorders
            linesplit = line.split()
            textid = int(re.search("\d+", linesplit[-1]).group(0))
            if textid == obid:
                x102 = float(linesplit[1].split(",")[0])
                y102 = float(linesplit[2].split(",")[0])
        reg102.close()
    if os.path.exists(path141_1reg) == 1:
        reg141 = open(path141_1reg, "r")
        x141, y141 = -1, -1
        for line in reg141:
            linesplit = line.split()
            textid = int(re.search("\d+", linesplit[-1]).group(0))
            if textid == obid:
                x141 = float(linesplit[1].split(",")[0])
                y141 = float(linesplit[2].split(",")[0])
        reg141.close()

    '''
    if load_image:
        if os.path.exists(path102) == 1:
            cmd = "xpaset -p ds9 frame 5"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path102
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path102_0reg
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path102_1reg
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to %f %f image" % (x102, y102)
            os.system(cmd)
        if os.path.exists(path141) == 1:
            cmd = "xpaset -p ds9 frame 6"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path141
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path141_0reg
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path141_1reg
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to %f %f image" % (x141, y141)
            os.system(cmd)
    '''
    # pan to the coordinates of this object
    if os.path.exists(path102):
        panDispersed(x102, y102, grism="G102")
    if os.path.exists(path141):
        panDispersed(x141, y141)


### AH: 1/26/28: I do not think we are using this?
def createAltGrismRegion(grism):
    workingdir = os.getcwd()
    par = os.path.dirname(workingdir)
    cat110 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F110.cat")
    cat140 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F140.cat")
    cat160 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F160.cat")

    if grism == "G102":
        if os.path.exists(cat110) == 1:
            cat = np.genfromtxt(
                cat110,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        else:
            print(cat110)
            return 0
    if grism == "G141":
        if os.path.exists(cat140) == 1:
            cat = np.genfromtxt(
                cat140,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        elif os.path.exists(cat160) == 1:
            cat = np.genfromtxt(
                cat160,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        else:
            print("nope2")
            return 0

    f = open(os.path.join(workingdir, grism + "_temp.reg"), "w")
    f.write(
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    )
    f.write("image\n")
    cenx = np.zeros(cat["num"].shape[0], dtype=float)
    ceny = np.zeros(cat["num"].shape[0], dtype=float)
    for i in range(cat["num"].shape[0]):
        objid = cat["num"][i]
        gfile = os.path.join(
            par, "%s_DRIZZLE/aXeWFC3_%s_mef_ID%i.fits" % (grism, grism, objid)
        )
        if os.path.isfile(gfile):
            hdr = fits.getheader(gfile, 1)
            # coords of bounding box
            boxx = np.array([hdr["bb0x"], hdr["bb1x"]])
            boxy = np.array([hdr["bb0y"], hdr["bb1y"]])
            cenx[i] = boxx[0] + (boxx[1] - boxx[0]) / 2.0
            ceny[i] = boxy[0] + (boxy[1] - boxy[0]) / 2.0
            slitwidth = hdr["slitwidt"]
            sx = 184
            sy = 8
            #            sy = slitwidth * cat['a_img'][i]
            f.write(
                "box(%f,%f,%i,%i,0) # text={%i}\n" % (cenx[i], ceny[i], sx, sy, objid)
            )
    f.close()
    return cenx, ceny


def panDirect(ra, dec, grism="F200W"):
    # Pan to coords in frame
    if grism == "F115W":
        fno = "4"
    elif grism == "F150W":
        fno = "5"
    elif grism == "F200W":
        fno = "6"
    cmd = "xpaset -p ds9 frame " + fno
    os.system(cmd)
    cmd = "xpaset -p ds9 pan to " + ra + " " + dec  # +' fk5'
    os.system(cmd)


def panDispersed(xx, yy, grism="G141"):  # MB
    # Pan to coords in frame
    if grism == "G141":
        fno = "6"
    else:
        fno = "5"
    cmd = "xpaset -p ds9 frame " + fno
    os.system(cmd)
    cmd = "xpaset -p ds9 pan to %f %f image" % (xx, yy)
    os.system(cmd)

## FH 2/11/25:
def reloadReg(path_to_data="",filter="F444W"):
    # workingdir = os.getcwd()
    # dirpts = workingdir.split("/")[1:-1]
    # par_root_dir = "/"
    # for pdir in dirpts:
    #     par_root_dir = par_root_dir + pdir + "/"

    ds9_title = "POPPIES_DIRECT"

    regfile = glob(path_to_data + "*_{}.reg".format(filter))
    # reload direct image region files
    if os.path.exists(regfile):
        cmd = f"xpaset -p ds9 {ds9_title} frame 1"
        os.system(cmd)
        cmd = (
            "xpaset -p ds9 regions file "
            + regfile
        )
        os.system(cmd)
    else:
        print("Region file not found for filter " + str(filter))
    # if os.path.exists(par_root_dir + "DATA/DIRECT_GRISM/F160_drz.reg") == 1:
    #     cmd = "xpaset -p ds9 frame 4"
    #     os.system(cmd)
    #     cmd = (
    #         "xpaset -p ds9 regions file "
    #         + par_root_dir
    #         + "DATA/DIRECT_GRISM/F160_drz.reg"
    #     )
    #     os.system(cmd)
    # elif os.path.exists(par_root_dir + "DATA/DIRECT_GRISM/F140_drz.reg") == 1:
    #     cmd = "xpaset -p ds9 frame 4"
    #     os.system(cmd)
    #     cmd = (
    #         "xpaset -p ds9 regions file "
    #         + par_root_dir
    #         + "DATA/DIRECT_GRISM/F140_drz.reg"
    #     )
    #     os.system(cmd)



## PASSAGE version:
def reloadReg_PASSAGE():
    workingdir = os.getcwd()
    dirpts = workingdir.split("/")[1:-1]
    par_root_dir = "/"
    for pdir in dirpts:
        par_root_dir = par_root_dir + pdir + "/"

    # reload direct image region files
    if os.path.exists(par_root_dir + "DATA/DIRECT_GRISM/F110_drz.reg") == 1:
        cmd = "xpaset -p ds9 frame 3"
        os.system(cmd)
        cmd = (
            "xpaset -p ds9 regions file "
            + par_root_dir
            + "DATA/DIRECT_GRISM/F110_drz.reg"
        )
        os.system(cmd)
    if os.path.exists(par_root_dir + "DATA/DIRECT_GRISM/F160_drz.reg") == 1:
        cmd = "xpaset -p ds9 frame 4"
        os.system(cmd)
        cmd = (
            "xpaset -p ds9 regions file "
            + par_root_dir
            + "DATA/DIRECT_GRISM/F160_drz.reg"
        )
        os.system(cmd)
    elif os.path.exists(par_root_dir + "DATA/DIRECT_GRISM/F140_drz.reg") == 1:
        cmd = "xpaset -p ds9 frame 4"
        os.system(cmd)
        cmd = (
            "xpaset -p ds9 regions file "
            + par_root_dir
            + "DATA/DIRECT_GRISM/F140_drz.reg"
        )
        os.system(cmd)
