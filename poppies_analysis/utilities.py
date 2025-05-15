#### This python script has a bunch of very nifty utility functions that are used ubiquitously

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import glob
import os
from astropy.table import Table
import astropy.io.ascii as asciitable
import pandas as pd


def gaussian(x, mu, sigma):
    return np.exp(-np.power(x - mu, 2.) / (2. * np.power(sigma, 2.))) # / (sigma * np.sqrt(2. * np.pi))



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# FH 2/4/25:
def read_config(config, availgrism='F444W'): 
    ''' availgrism is the filter(s) '''
    configfile = open(config, 'r')  
    config_pars = {} 
    for line in configfile:
        if ( (line[0] != '#') & (len(line)>1)): 
            name = line.split()[0] 
            
            if name == 'node_wave': 
                tmpnode = [float(s) for s in line.split()[1::]]
                nodelam = [x for x in tmpnode if 
                ((x >= config_pars['lambda_min_{}'.format(availgrism)]) and 
                (x <= config_pars['lambda_max_{}'.format(availgrism)]))]

                # if availgrism == 'F277W':
                #     nodelam = [x for x in tmpnode if x <= config_pars['transition_wave']]
                # elif availgrism.lower() == 'g141':
                #     nodelam = [x for x in tmpnode if x > config_pars['transition_wave']]
                # else:
                #     nodelam = tmpnode

                config_pars.setdefault('node_wave', []) 
                for l in nodelam:  
                    config_pars['node_wave'].append(l)

            elif name == 'mask_region1':

                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region1', []) 
                config_pars['mask_region1'].append(masklam[0]) 
                config_pars['mask_region1'].append(masklam[1]) 
            elif name == 'mask_region2': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region2', []) 
                config_pars['mask_region2'].append(masklam[0]) 
                config_pars['mask_region2'].append(masklam[1])
            elif name == 'mask_region3': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region3', []) 
                config_pars['mask_region3'].append(masklam[0]) 
                config_pars['mask_region3'].append(masklam[1])  
            elif name == 'mask_region4': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region4', []) 
                config_pars['mask_region4'].append(masklam[0]) 
                config_pars['mask_region4'].append(masklam[1])  
            elif name == 'mask_region5': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region5', []) 
                config_pars['mask_region5'].append(masklam[0]) 
                config_pars['mask_region5'].append(masklam[1])  
            elif name == 'mask_region6': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region6', []) 
                config_pars['mask_region6'].append(masklam[0]) 
                config_pars['mask_region6'].append(masklam[1])  
            elif name == 'mask_region7': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region7', []) 
                config_pars['mask_region7'].append(masklam[0]) 
                config_pars['mask_region7'].append(masklam[1])  
            elif name == 'mask_region8': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region8', []) 
                config_pars['mask_region8'].append(masklam[0]) 
                config_pars['mask_region8'].append(masklam[1])  
            else:  
                val  = line.split()[1] 
                if is_number(val) == True:  val = float(val)
                if val == 'True': val = True
                if val == 'False': val = False 
                config_pars[line.split()[0]] = val 
    configfile.close()
    return config_pars    



# Creates .reg files for direct and dispersed images
# Written by Mason Huberty and adapted by KVN to work as part of the line finding software
# modified by FH for POPPIES NIRCam data

def write_obj_region(parno, path_to_data, catalog, regfile_name, xoffset, yoffset, w, b_width, b_length): 
    file = open(path_to_data + "/POPPIES" + str(parno) + "/DATA/" + "POPPIES" + str(parno) + regfile_name, 'a')
    for i in range(len(catalog)):
        ra, dec = catalog['RA'][i], catalog['DEC'][i]
        x, y = w.all_world2pix(ra, dec, 1)
        file.write("box("+str(x-xoffset)+','+str(y-yoffset)+','+str(b_width)+','+str(b_length)+',0.0) # color=red text={'+str(catalog['NUMBER'][i])+'} font="times 10 bold italic" textangle=30\n')
    file.close


# NOTE - This is for NIRCam  
#        F.H. started 12/17/24
def create_regions(parno, path_to_data, filters):

    # spec_cat =  glob.glob(path_to_data + "/POPPIES" + str(parno) + "/*_A_*_i2d.cat") 
    
    # all_objects = []

    ## For now, just making separate region files for separate filters
    for filt in filters:
            
        phot_cat = glob.glob(path_to_data + "POPPIES" + str(parno) + "/*_{}_i2d.cat".format(str(filt))) 

        cat=asciitable.read(phot_cat[0])

        # Direct image region files
        f = open(path_to_data + "/POPPIES" + str(parno) + "/POPPIES" + str(parno) + 'regions_phot_{}.reg'.format(str(filt)),'a')
        for i in range(len(cat)):
            f.write("WCS;circle("+str(cat['RA'][i])+','+str(cat['DEC'][i])+',0.5") # color=green text={'+str(cat['NUMBER'][i])+'} font="times 10 bold italic" textangle=30\n')
        f.close()

    #This and subsequent are for the first order beams. Offsets taken from the config files
    # KVN :  Adding code to first check if the paths exist. Only create region file if filter/orientation is available 
    # FH updated 12/20/24

    ## for ref: write_obj_region(parno, path_to_data, catalog, regfile_name, xoffset, yoffset, w, b_width, b_length)
    
    # f277grism_R = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/spec2D/*_F277W_*_i2d.fits")
    # if len(f277grism_R) != 0:
    #     write_obj_region(parno, path_to_data, cat, "F277r_grism.reg", 0.6548681566074263, 33.73739138173772, w = WCS(f277grism_R[0][1]), 
    #                     b_width = 10.0, b_length = 93.54)
    
    # f115grism_C = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f115w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
    # if len(f115grism_C) != 0:
    #     write_obj_region(parno, path_to_data, cat, "F115c_grism.reg", 31.91107156101387, 1.3922939626209256, w = WCS(f115grism_C[0]), 
    #                     b_width = 97.28751330105166, b_length = 10.0)

    # f150grism_R = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f150w*gr150r_drz_sci.fits") # Added KVN 19-Aug-2024
    # if len(f150grism_R) != 0:
    #     write_obj_region(parno, path_to_data, cat, "F150r_grism.reg", 0.6548681566074263, 106.79254657227568, w = WCS(f150grism_R[0]), 
    #                     b_width = 10.0, b_length = 93.54)
    
    # f150grism_C = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f150w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
    # if len(f150grism_C) != 0:
    #     write_obj_region(parno, path_to_data, cat, "F150c_grism.reg", 96.44444, 0.6548681566074263, w = WCS(f150grism_C[0]), 
    #                     b_width = 93.54, b_length = 10.0)
    

    ### FH started this 2/3/25 but leaving for now:

    # if filt == "F444W":

    #     f444grism = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DirectImages/*_F444W_i2d.fits") # Added FH 12/23/24
    #     if len(f444grism) != 0:    
    #         write_obj_region(parno, path_to_data, cat, "F444r_grism.reg", 765.438247, 9.66479920022, w = WCS(f444grism[0][1]), 
    #                         b_width = 10.0, b_length = 131.78)        
    
    #     if len(f444grism) != 0:
    #         write_obj_region(parno, path_to_data, cat, "F444c_grism.reg", 24.1334945641, 765.438247, w = WCS(f444grism[0][1]), 
    #                         b_width = 127.806, b_length = 10.0)        

    # elif filt == "F277W":

    #     f444grism = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DirectImages/*_F277W_i2d.fits") # Added FH 12/23/24
    #     if len(f444grism) != 0:    
    #         write_obj_region(parno, path_to_data, cat, "F444r_grism.reg", 765.438247, 9.66479920022, w = WCS(f444grism[0][1]), 
    #                         b_width = 10.0, b_length = 131.78)        
    
    #     if len(f444grism) != 0:
    #         write_obj_region(parno, path_to_data, cat, "F444c_grism.reg", 24.1334945641, 765.438247, w = WCS(f444grism[0][1]), 
    #                         b_width = 127.806, b_length = 10.0)    


def add_header_keyword(parno, path_to_data):

    main_directory = os.path.join(path_to_data, f"Par{parno:s}")
    files = glob.glob(os.path.join(main_directory, 'spec2D/*.fits'))
    test_file = fits.open(files[0])
    header = test_file[2].header

    if header['RADESYS'] != 'ICRS':
        print('Headers are not updated. Updating all now.')
        for j in files:
            # Open the file header for viewing and load the header
            hdulist = fits.open(j)
            for i in range(len(hdulist)):
                header = hdulist[i].header
                try: header['RADESYS'] = 'ICRS'
                except: print('no RADESYS in header')

            hdulist.writeto(j, overwrite='True')


## adapting for POPPIES, FH 1/23/25
def make_spectra_dat_files(parno, path_to_data, create_files_anyway = False):

    main_directory = os.path.join(path_to_data, f"POPPIES{parno:s}")
    spec_directory = os.path.join(main_directory, "Spec1D2D")

    # os.system(f"mkdir -p {os.path.join(main_directory, 'Spectra'):s}")

    files = sorted(glob.glob(os.path.join(spec_directory, '*.fits')))


    # Check if converted files already exist. If they do not, carry on with the conversion.
    # Otherwise, this step can be skipped.
    # Here, I just check if there are more converted files than objects (there can be 1-3 per object)
    # depending on how many filters are available for each object
    already_converted_files = glob.glob(main_directory + '/Spectra/*.dat')

    # print(len(files),len(already_converted_files))

    if (len(files) >= len(already_converted_files)) or create_files_anyway == True:
        # print(len(already_converted_files))
        # print(len(files))

        print('\nThere are ' + str(len(files) - len(already_converted_files)) + ' POPPIES files to convert.\n')


    for f in files:

        if f not in already_converted_files:
                
            try:
                fff = fits.open(f,ignore_missing_simple=True)

                # for ext in range(1, len(fff)):

                ext = 1 ## coresponds to 1D spec for a give filter

                tb = Table.read(fff[ext], format='fits').to_pandas()
                t_out = pd.DataFrame({})

                ### !!! IMPORTANT !!!
                ### VM: This is temporary while we fix this in the pipeline
                ### The R/C spectra have already been treated for the following operations
                # if "EXTVER" not in fff[ext].header:
                    
                t_out['wave'] = tb['WAVELENGTH'] * 10000 #micron to A
                
                # t_out['flux'] = tb['FLUX']
                # t_out['error'] = tb['FLUXERR']

                t_out['flux'] = tb['OPTFLUX']
                t_out['error'] = tb['OPTFLUXERR']

                t_out['contam'] = tb['CONTAMFLUX']
                t_out['zeroth'] = np.zeros(len(tb['WAVELENGTH'])).astype('int')
                
                # t_out = Table.from_pandas(t_out)
                # t_out = t_out.filled(0.0) # Replace nans with zeros
                # # Spectra dispersed beyond the chip have zero fluxes that must be replaced to prevent crashes in fitting.
                # t_out['flux'][np.where(t_out['flux'] == 0.0)] = np.median(t_out['flux'][np.where(t_out['flux'] != 0.0)])
                # t_out['error'][np.where(t_out['error'] == 0.0)]=np.median(t_out['error'][np.where(t_out['error'] != 0.0)])

                t_out = Table.from_pandas(t_out)
                ### VM: Explicitly fix nans to 0s
                for col in t_out.columns:
                    t_out[col][np.isnan(t_out[col])] = 0
                t_out = t_out.filled(0.0) # Replace nans with zeros

                # Spectra dispersed beyond the chip have zero fluxes that must be replaced to prevent crashes in fitting.
                t_out['flux'][np.where(t_out['flux'] == 0.0)] = np.median(t_out['flux'][np.where(t_out['flux'] != 0.0)])
                t_out['error'][np.where(t_out['error'] == 0.0)] = np.median(t_out['error'][np.where(t_out['error'] != 0.0)])

                # t_out['flux'][np.where(np.isnan(t_out['flux']))] = np.median(t_out['flux'][np.where(t_out['flux'] != 0.0)])
                # t_out['flux'][np.where(np.isnan(t_out['flux']))] = np.median(t_out['flux'][np.where(t_out['error'] != 0.0)])

                # FH 2/5/25 - catch those that have wavelength = 0 at the end:
                if t_out['wave'][-1] == 0:
                    t_out.remove_row(-1)

                print('saving ' + str(f) + ' to '+ str(main_directory) + '/Spectra/')

                # for filt in ["115", "150", "200"]:
                        
                t_out.write(main_directory + '/Spectra/'+os.path.basename(f).replace('.fits', '_1D.dat'), 
                        format='ascii.fixed_width_two_line', overwrite=True) 

            except Exception as e:
                print('Could not process {}, due to error {}'.format(f,e))
                print(t_out.info())
                pass


## FH 2/3/25
def find_filters(path_to_data,par):
    ''' searches catalogs for list of filters - and orientations '''
    
    cats = glob.glob(path_to_data + "POPPIES" + str(par) + "/*_i2d.cat")

    filts = []

    for i in cats:
        print(i)
        y = str(i).split('_i2d.cat')[0].split('_F')[-1]
        filt = str('F')+str(y)
        filts.append(filt)

    unique_filts = np.unique(filts)

    outfile = open(path_to_data+"POPPIES" + str(par)+'/POPPIES{}_filters.dat'.format(par), 'w')

    outfile.write('filter ' + 'orientation' + '\n')

    for i in unique_filts:

        Rfiles = glob.glob(path_to_data + "POPPIES" + str(par) + "/Spec1D2D/*_{}_R_*.fits".format(i))

        Cfiles = glob.glob(path_to_data + "POPPIES" + str(par) + "/Spec1D2D/*_{}_C_*.fits".format(i))

        if ((len(Rfiles) != 0) and (len(Cfiles) != 0)):

            outfile.write(str(i) + '  ' + 'R+C' + '\n')

        elif ((len(Rfiles) != 0) and (len(Cfiles) == 0)):

            outfile.write(str(i) + '  ' + 'R' + '\n')

        elif ((len(Rfiles) == 0) and (len(Cfiles) != 0)):

            outfile.write(str(i) + '  ' + 'C' + '\n')


## FH 1/30/25
def make_full_list(path_to_data,path_to_out,par,filters,verbose=True):
    ''' makes full object list from se cats'''

    #check that linelist directory exists in Output:
    if not os.path.exists(path_to_out+'/linelist/'):
        
        command = "mkdir {} ".format(str(path_to_out+'/linelist/'))
        os.system(command)

    # file for full object list:
    outfile = open(path_to_out+'/linelist/POPPIES'+str(par) + 'objects.dat', 'w')

    # find all the catalogs:
    for filt in filters:

        secats = glob.glob(path_to_data + "POPPIES" + str(par) + "/*_{}_i2d.cat".format(str(filt))) 
        
        secats.sort()

        for j in secats:

            y = str(j).split('_i2d.cat')[0].split('_F')[-1]
            filt = str('F')+str(y)

            cat = asciitable.read(j)

            if verbose == True:
                print("I found the following photometric catalogs...\n") 
                print(j)  

            ids = cat['NUMBER']

            for id in ids:
                # print(filt,id)
                outfile.write(str(f'{par:03}') + '  ' + str(filt) + '  ' + str(id) + '\n')

        
    outfile.close()
        # ids_all.append(ids)

    # ids_unique = np.unique(ids)

    # ids_b = cat_b['NUMBER']
    
    # filts = [str(filt) for i in ids]

    # # file for full object list:
    # outfile = open(path_to_out+'/linelist/Par'+str(par) + 'objects_{}.dat'.format(str(filt)), 'w')

    # for id in ids:
    #     outfile.write(str(par) + '  ' + str(filt) + '  ' + str(id) + '\n')


## FH 2/3/25
def make_file_structure(path_to_data,par):

    spec1d2dpath = path_to_data + "/POPPIES" + str(par) + "/Spec1D2D/"
    spectrapath = path_to_data + "/POPPIES" + str(par) + "/Spectra/"
    directpath = path_to_data + "/POPPIES" + str(par) + "/DirectImages/"

    spec1d2d = glob.glob(spec1d2dpath)
    
    spectra = glob.glob(spectrapath) 

    direct = glob.glob(directpath)

    if len(spec1d2d) == 0:
        print('Creating Spec1D2D directory')
        os.mkdir(spec1d2dpath)

    if len(spectra) == 0:
        print('Creating Spectra directory')
        os.mkdir(spectrapath)

    if len(direct) == 0:
        print('Creating DirectImage directory')
        os.mkdir(directpath)

    spec1d2dfiles = glob.glob(path_to_data + "/POPPIES" + str(par) + "/jw*_ext/*.fits")

    if len(spec1d2dfiles) != 0:
        print('Moving 1D and 2D spec files to correct path')

        for file in spec1d2dfiles:

            command = "mv {} {} ".format(file,spec1d2dpath)

            os.system(command)

    directfiles = glob.glob(path_to_data + "/POPPIES" + str(par) + "/*_i2d.fits")


    if len(directfiles) != 0:
        print('Moving direct image files to correct path')

        for file in directfiles:

            command = "mv {} {} ".format(file,directpath)

            os.system(command)

## FH 2/4/25:
def quick_flux_max(wave,flux,err,wavemin,wavemax):
    
    ''' very simple function to find maximum flux of a spectrum '''
    # config_pars = read_config(config, filter)

    # wavemin,wavemax = config_pars['lambda_min_{}'.format(filter)], config_pars['lambda_max_{}'.format(filter)]
    minwave_offset = 300  #small offset
    sn = flux/err #S/N array

    snsorted = np.sort(sn) #sort S/N array
    
    id = -1
    maxind = np.where(sn == snsorted[id])[0][0]  #Max S/N index

    #iterate towards smaller S/N if needed until we reach a wavelength that is within filter limits
    while (wave[maxind] < wavemin+minwave_offset) or (wave[maxind] > wavemax-minwave_offset):
        id += -1
        try:
            maxind = np.where(sn == snsorted[id])[0][0]

        except Exception as e: #fail-safe to prevent crashing
            print("Max. S/N didn't work: \n", e)
            maxind = 0


    wave_maxind = wave[maxind]
    sn_maxind = sn[maxind]

    return wave_maxind, sn_maxind



# # NOTE - The x & y offset values are specific for NIRISS. 
# #        These will need to be updated for NIRCam data!
# def create_regions_OLD(parno, path_to_data):

#     spec_cat = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/DIRECT_GRISM/Par*spec*.fits")
#     hdul = fits.open(spec_cat[0])
#     cat=hdul[1].data

#     # Direct image region files
#     f = open(path_to_data + "/POPPIES" + str(parno) + "/DATA/" + "POPPIES" + str(parno) + 'regions_phot.reg','a')
#     for i in range(len(cat)):
#         f.write("WCS;circle("+str(cat['ra'][i])+','+str(cat['dec'][i])+',0.5") # color=green text={'+str(cat['id'][i])+' z='+str(round(cat['redshift'][i],3))+'} font="times 10 bold italic" textangle=30\n')
#     f.close()

#     #This and subsequent are for the first order beams. Offsets taken from the config files
#     # KVN :  Adding code to first check if the paths exist. Only create region file if filter/orientation is available 
#     f115grism_R = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f115w*gr150r_drz_sci.fits")
#     if len(f115grism_R) != 0:
#         write_obj_region(parno, path_to_data, cat, "F115r_grism.reg", 0.6548681566074263, 33.73739138173772, w = WCS(f115grism_R[0]), 
#                         b_width = 10.0, b_length = 93.54)
    
#     f115grism_C = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f115w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
#     if len(f115grism_C) != 0:
#         write_obj_region(parno, path_to_data, cat, "F115c_grism.reg", 31.91107156101387, 1.3922939626209256, w = WCS(f115grism_C[0]), 
#                         b_width = 97.28751330105166, b_length = 10.0)

#     f150grism_R = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f150w*gr150r_drz_sci.fits") # Added KVN 19-Aug-2024
#     if len(f150grism_R) != 0:
#         write_obj_region(parno, path_to_data, cat, "F150r_grism.reg", 0.6548681566074263, 106.79254657227568, w = WCS(f150grism_R[0]), 
#                         b_width = 10.0, b_length = 93.54)
    
#     f150grism_C = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f150w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
#     if len(f150grism_C) != 0:
#         write_obj_region(parno, path_to_data, cat, "F150c_grism.reg", 96.44444, 0.6548681566074263, w = WCS(f150grism_C[0]), 
#                         b_width = 93.54, b_length = 10.0)
    
#     f200grism_R = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f200w*gr150r_drz_sci.fits") # Added KVN 19-Aug-2024
#     if len(f200grism_R) != 0:    
#         write_obj_region(parno, path_to_data, cat, "F200r_grism.reg", 0.6548681566074263, 204.8370874255101, w = WCS(f200grism_R[0]), 
#                         b_width = 10.0, b_length = 131.78)        
    
#     f200grism_C = glob.glob(path_to_data + "/POPPIES" + str(parno) + "/DATA/*f200w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
#     if len(f200grism_C) != 0:
#         write_obj_region(parno, path_to_data, cat, "F200c_grism.reg", 200.9228, 0.6548681566074263, w = WCS(f200grism_C[0]), 
#                         b_width = 127.806, b_length = 10.0)        
    

