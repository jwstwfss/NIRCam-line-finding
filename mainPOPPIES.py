##### This is the main python file for NIRCAM-line-finding #####
### This version started 1/30/25
## Farhan Hasan (fhasan@stsci.edu)

import poppies_analysis as poppies
from poppies_analysis import utilities
import os
import re
import glob
import astropy.io.ascii as asciitable
import time


# # Open two ds9 windows:
# os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title POPPIES_DIRECT &')
# os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title POPPIES_spec2D &')


### Please update these directories to match yours:
CODE_DIR = "/Users/fhasan/Desktop/Research_STScI/POPPIES/line-finding/poppies_analysis"
OUTPUT_DIR = "/Users/fhasan/Desktop/Research_STScI/POPPIES/Output"
DATA_DIR = "/Users/fhasan/Desktop/Research_STScI/POPPIES/Data/"


#############################################################
####### You should not need to change anything below. #######
#############################################################

print("Code started at: ", time.gmtime())

parno = input('\033[94m' + "Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
# pull out only the number, in case the user entered e.g. "Par1"
while True:
    try:
        parno = re.findall(r'\d+', str(parno))[0]
    except: 
        parno = input('\033[94m' + "A parallel field number is required. Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
        continue
    else:
        break


if __name__ == "__main__":
    
    ## FH 2/3/25: make file structure if one does not exist for the data:
    utilities.make_file_structure(DATA_DIR,str(parno))
   
    # FH 2/3/25: check for how many filters there are for this field:
    filterfile = glob.glob(DATA_DIR + "Par" + str(parno) + "/Par{}_filters.dat".format(parno))
   
    # FH 2/3/25: make list of filters if one does not exist
    if len(filterfile) == 0:
        utilities.find_filters(DATA_DIR,str(parno))

        filterfile = glob.glob(DATA_DIR + "Par" + str(parno) + "/Par{}_filters.dat".format(parno))

    filterlist = asciitable.read(filterfile[0])


    #FH 1/30/25: check for full object list from SE cats
    objectfiles = glob.glob(OUTPUT_DIR + '/linelist/Par'+str(parno) + 'objects.dat')
    
    if len(objectfiles) == 0:

    # for filt in filterlist['filter']:
        
        utilities.make_full_list(DATA_DIR,OUTPUT_DIR,str(parno),filterlist['filter'])

        objectfiles = glob.glob(OUTPUT_DIR + '/linelist/Par'+str(parno) + 'objects.dat')

    objectlist = asciitable.read(objectfiles[0],names=['parno','filter','id'])

    # Check if .dat spec files exists or not
    specdatfiles = glob.glob(DATA_DIR + "Par" + str(parno) + "/Spectra/*.dat")
    
    if len(specdatfiles) == 0:
        # utilities.add_header_keyword(parno = parno, path_to_data = DATA_DIR) ## this one is commented out until RADESYS headers are provided in the reduced data

        try:
            utilities.make_spectra_dat_files(parno = parno, path_to_data = DATA_DIR)

        except Exception as e:
            print('Error making spec dat files due to error: ', e )
            pass


    # check if region files exist. If not, run code to create necessary region files
    regionfiles = glob.glob(DATA_DIR + "/Par" + str(parno) + "/*.reg")
    
    if len(regionfiles) == 0:
        print('\033[94m' + "No region files found, creating those for you now."  + '\033[0m')
        utilities.create_regions(parno, DATA_DIR, filterlist['filter'])


    # move to the directory where you want your outputs.
    os.chdir(OUTPUT_DIR)

    ## FH 2/3/25: ask user if they prefer CWT or full object list:

    # objtypes = input("Please hit ENTER to go through all objects identified by SE, or 'c' for wavelength decomposition \n> ")

    objtypes = input("Please hit ENTER or 'c' to perform wavelength decomposition, or 'all' to go through all objects identified by SE \n> ")

    if objtypes.strip().lower() == "all":

        linelist = objectlist
    
        #run the measure_z_interactive code to go through objects individually:
        poppies.measure_z_interactive_all(path_to_data=DATA_DIR, path_to_code=CODE_DIR, parno=parno)


    elif (objtypes == "") or (objtypes.strip().lower() == "c"):

        # check if line list exists. If not, run code to create the linelist
        print(OUTPUT_DIR + "/linelist/Par"+str(parno)+"lines.dat")
        linelist = glob.glob(OUTPUT_DIR + "/linelist/Par"+str(parno)+"lines.dat")
        
        if len(linelist) == 0:
            print('\033[94m' + "No line list file found, creating the line list for you now." + '\033[0m')
            
            poppies.loop_field_cwt(path_to_data=DATA_DIR, path_to_code=CODE_DIR, parno=parno)

        #run the measure_z_interactive code to go through objects individually:
        poppies.measure_z_interactive(path_to_data=DATA_DIR, path_to_code=CODE_DIR, parno=parno)
       


