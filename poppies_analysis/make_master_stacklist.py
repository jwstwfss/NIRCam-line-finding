#### this code makes subsamples for using with the stacking code. a host of various cuts will be included.
### run in directory where the catalogs live.


import astropy.io.fits as fits
import astropy.io.ascii as asciitable
import numpy as np
import os
import copy
#from __future__ import print_function #M.D.R#


########### STEP ZERO ###### CATALOGS

####### WISP CATALOGS AND PATHS
### read files
#M.D.R# masscat_wisp = asciitable.read('mederr_MARC_avgal_newparam.dat')
#masscat_wisp = asciitable.read('MUDF_SPS_20210416.txt') #M.D.R#
masscat_wisp = asciitable.read('MUDF_SPS_20210416_manual_hst.txt') #M.D.R# Merge of 'MUDF_SPS_20210416.txt' and 'MUDF_SPS_HSTsources.cat'.
# The file above uses a place holder redshift of 1.50 for all grism-only sources as I neglected to include these in my list of IDs to Matteo.
#M.D.R# hdu3 = fits.open('selected_photometry_catalog.fits')
#M.D.R# path_wisp = '/Volumes/Thunderbay/wisps/mzr_refit/'
path_wisp = '/Users/mrevalski/Documents/stsci/research/mudf/analysis/metallicity/run07/wisp_analysis'

## extra steps for wisp cats.
#M.D.R# wispcat = hdu3[1].data   ### this extra list is necessary because the par numbers are not in masscat_wisp


####### 3DHST CATALOGS AND PATHS
#M.D.R# hdu1 = fits.open('selected_c_cami.fits')
#M.D.R# path_3dhst = '/Volumes/Thunderbay/3DHST/mzr_refit/'


### agn matches
#M.D.R# xray = asciitable.read('../agn/xray_matches', format = 'no_header')
#M.D.R# donley = asciitable.read('../agn/donley_agn', format = 'no_header')


#### extra steps for 3dhst cats.
#M.D.R# masscat_3dhst = hdu1[1].data
#M.D.R# w=np.where(masscat_3dhst['field'] == 'GOODS-S')
#M.D.R# masscat_3dhst['field'][w] = 'GOODSS'
#M.D.R# w=np.where(masscat_3dhst['field'] == 'GOODS-N')
#M.D.R# masscat_3dhst['field'][w] = 'GOODSN'
#M.D.R# fieldnames_cat = masscat_3dhst['field']
#M.D.R# objid = masscat_3dhst['phot_id']
####  we're going to have to do this one field at a time.
### gather the fields, preserve them in the order in which they appear.
#M.D.R# dummy = np.unique(masscat_3dhst['field'], return_index = True)
#M.D.R# dummy2 = np.sort(dummy[1])
#M.D.R# fields = masscat_3dhst['field'][dummy2]


##### STEP ONE --- remove objects from these catalogs that were removed from the sample at some stage, or rejected in the fitting.
##############################
##############################
#### 3dhst ####################
##############################
##############################

#M.D.R# keepers = []
#M.D.R# z_3dhst  = []
#M.D.R# foiii_3dhst = []
#M.D.R# eoiii_3dhst = []
#M.D.R# fhb_3dhst = []
#M.D.R# ehb_3dhst = []
#M.D.R# fha_3dhst = []
#M.D.R# eha_3dhst = []
#M.D.R# fsii_3dhst = []
#M.D.R# esii_3dhst = []
#M.D.R# foii_3dhst = []
#M.D.R# eoii_3dhst = []
#M.D.R# ew_oiii_obs_3dhst  = []


#M.D.R# for fieldname in fields:
    #### first define the catalog and comments files:
#M.D.R#     if os.path.exists(path_3dhst+  '/' + fieldname + '_output_alaina-mzr/'):
#M.D.R#         catalog  = path_3dhst + '/' + fieldname + '_output_alaina-mzr/' + fieldname + '_catalog_alaina-mzr.dat'
#M.D.R#         comments = path_3dhst + '/' + fieldname + '_output_alaina-mzr/' + fieldname + '_comments_alaina-mzr.dat'
#M.D.R#     elif os.path.exists(path_3dhst + '/' + fieldname + '_output_marc-mzr/'):
#M.D.R#         catalog  = path_3dhst + '/' + fieldname + '_output_marc-mzr/' + fieldname + '_catalog_marc-mzr.dat'
#M.D.R#         comments = path_3dhst + '/' + fieldname + '_output_marc-mzr/' + fieldname + '_comments_marc-mzr.dat'
#M.D.R#     else:
#M.D.R#         comments = None
#M.D.R#         catalog = None
#M.D.R#         print 'Could not find fit data directory for ' + fieldname

    #### now read the catalog and the comments files for this field
#M.D.R#     catalog_data =asciitable.read(catalog, format = 'sextractor')
#M.D.R#     catalog_objid = catalog_data['ObjID']

    ### the comments file is not easily readable.
#M.D.R#     f = open(comments)
#M.D.R#     obj_comments = []
#M.D.R#     reject_list = []

#M.D.R#     for line in f:
#M.D.R#        x = line.split()
#M.D.R#        obj_comments.append(x[1])

#M.D.R#        if len(x) > 2 :
#M.D.R#           if x[2] == 'rejected':
#M.D.R#                 reject_list.append('rejected')
#M.D.R#           else:
#M.D.R#                 reject_list.append('keep')
#M.D.R#        else:
#M.D.R#             reject_list.append('keep')
#M.D.R#     f.close()


#M.D.R#     ind_field= np.where(fieldnames_cat == fieldname)
#M.D.R#     ind_field = ind_field[0]
#M.D.R#     obj_field = objid[ind_field]


#M.D.R#     for i in np.arange(len(obj_field)) :
#M.D.R#         if os.path.exists(path_3dhst+  '/' + fieldname + '_output_alaina-mzr/'):
#M.D.R#             specfile = path_3dhst + '/' + fieldname + '_output_alaina-mzr/fitdata/' + fieldname+  '_' +  '{:05d}'.format(int(obj_field[i])) + '_fitspec.dat'
#M.D.R#         elif os.path.exists(path_3dhst + '/' + fieldname + '_output_marc-mzr/'):
#M.D.R#             specfile = path_3dhst + '/' + fieldname + '_output_marc-mzr/fitdata/' + fieldname + '_' + '{:05d}'.format(int(obj_field[i])) + '_fitspec.dat'
#M.D.R#         else:
#M.D.R#             specfile = '~/path/to/nonsense/'


#M.D.R#         if os.path.exists(specfile) == True :
#M.D.R#             check1 = True
#M.D.R#         else :
#M.D.R#             check1 = False

#M.D.R#         if (int(obj_field[i]) in catalog_objid) and ( '{:05d}'.format(int(obj_field[i]))  in obj_comments) :
#M.D.R#             check2 = True
#M.D.R#         else :
#M.D.R#             check2 = False


#M.D.R#         obj_comments = np.array(obj_comments)
#M.D.R#         w=np.where(obj_comments == '{:05d}'.format(int(obj_field[i]))   )

#M.D.R#         if np.size(w[0]) == 0 :
#M.D.R#             check3 = False
#M.D.R#         else:
#M.D.R#             if reject_list[w[0][0]] == 'keep' :
#M.D.R#                 check3 = True
#M.D.R#             else:
#M.D.R#                 check3 = False

        #print fieldname, '{:05d}'.format(int(obj_field[i])), check1, check2, check3

#M.D.R#         if (check1 == True) & (check2 == True) & (check3 == True) :
#M.D.R#             keepers.append(1)
#M.D.R#             w=np.where(catalog_objid == int(obj_field[i]))
#M.D.R#             z_3dhst.append(catalog_data['redshift'][w[0][-1]])  ###  just in case the object is in the catlag more than once, take the last one.
#M.D.R#             foiii_3dhst.append(catalog_data['oiii_flux'][w[0][-1]])
#M.D.R#             eoiii_3dhst.append(catalog_data['oiii_err'][w[0][-1]])
#M.D.R#             fha_3dhst.append(catalog_data['hanii_flux'][w[0][-1]])
#M.D.R#             eha_3dhst.append(catalog_data['hanii_err'][w[0][-1]])
#M.D.R#             fhb_3dhst.append(catalog_data['hb_flux'][w[0][-1]])
#M.D.R#             ehb_3dhst.append(catalog_data['hb_err'][w[0][-1]])
#M.D.R#             fsii_3dhst.append(catalog_data['sii_flux'][w[0][-1]])
#M.D.R#             esii_3dhst.append(catalog_data['sii_err'][w[0][-1]])
#M.D.R#             foii_3dhst.append(catalog_data['oii_flux'][w[0][-1]])
#M.D.R#             eoii_3dhst.append(catalog_data['oii_err'][w[0][-1]])
#M.D.R#             ew_oiii_obs_3dhst.append(catalog_data['oiii_EW_obs'][w[0][-1]])



#M.D.R#         else:
#M.D.R#             keepers.append(0)



#### LASTLY, cull the parent catalog.

#M.D.R# print 'Size of Parent 3dhst catalog: '
#M.D.R# print np.size(masscat_3dhst['field'])

#M.D.R# print 'Size of retained 3dhst catalog'
#M.D.R# print np.sum(keepers)


#M.D.R# keepers = np.array(keepers)
#M.D.R# w=np.where(keepers == 1)
#### gather what we need from these catalogs:
#M.D.R# id_3dhst = masscat_3dhst['phot_id'][w]
#M.D.R# field_3dhst = masscat_3dhst['field'][w]
#M.D.R# logm_3dhst = masscat_3dhst['logM_50_cami'][w]


#M.D.R# z_3dhst = np.array(z_3dhst)
#M.D.R# foiii_3dhst = np.array(foiii_3dhst)
#M.D.R# eoiii_3dhst = np.array(eoiii_3dhst)

#M.D.R# print np.size(z_3dhst)
#M.D.R# print np.size(id_3dhst)

##############################
##############################
#### WISP ####################
##############################
##############################
objid_parent = masscat_wisp['ID']   ####  integer numpy array
#M.D.R# par_parent = wispcat['par']
par_parent = [0] #M.D.R#


keepers = []
z_wisp = []
foiii_wisp = []
eoiii_wisp = []
fha_wisp = []
eha_wisp = []
fhb_wisp = []
ehb_wisp = []
fsii_wisp = []
esii_wisp = []
foii_wisp = []
eoii_wisp  = []

ew_oiii_obs_wisp = []





#M.D.R# for i  in np.arange(len(par_parent)):
for i  in np.arange(len(objid_parent)): #M.D.R#

#M.D.R#     if os.path.exists(path_wisp +  '/Par' + str(par_parent[i]) + '_output_alaina-mzr/') :
#M.D.R#             specfile = path_wisp + '/Par' + str(par_parent[i]) + '_output_alaina-mzr/fitdata/Par' + str(par_parent[i]) + '_BEAM_' + str(objid_parent[i]) + '_fitspec.dat'
#M.D.R#             catalog  = path_wisp  + '/Par' + str(par_parent[i])+ '_output_alaina-mzr/Par' + str(par_parent[i]) + 'list_catalog_alaina-mzr.dat'
#M.D.R#             comments = path_wisp  + '/Par' + str(par_parent[i])+ '_output_alaina-mzr/Par' + str(par_parent[i]) + 'list_comments_alaina-mzr.dat'

#M.D.R#     elif  os.path.exists(path_wisp +  '/Par' + str(par_parent[i]) + '_output_marc-mzr/') :
#M.D.R#             specfile = path_wisp + '/Par' + str(par_parent[i]) + '_output_marc-mzr/fitdata/Par' + str(par_parent[i]) + '_BEAM_' + str(objid_parent[i]) + '_fitspec.dat'
#M.D.R#             catalog  = path_wisp  + '/Par' + str(par_parent[i])+ '_output_marc-mzr/Par' + str(par_parent[i]) + 'list_catalog_marc-mzr.dat'
#M.D.R#             comments = path_wisp  + '/Par' + str(par_parent[i])+ '_output_marc-mzr/Par' + str(par_parent[i]) + 'list_comments_marc-mzr.dat'

    if os.path.exists(path_wisp +  '/Par0/Spectra/Par0_output_run07/') :
            specfile = path_wisp + '/Par0/Spectra/Par0_output_run07/' + '/fitdata/' + 'Par0_' + str(objid_parent[i]) + '_fitspec.dat'
            catalog = path_wisp  + '/Par0/Spectra/Par0_output_run07/' + 'mudf_line_flux_catalog_01.dat'
            comments = path_wisp  + '/Par0/Spectra/Par0_output_run07/' + 'mudf_line_flux_comments_01.dat'

    else:

        specfile = '~/path/to/nonsense/'


    ### check 1, does the spectrum exist
    if os.path.exists(specfile) == True :
            check1 = True
            #print('\ncheck1 passed!\n') #M.D.R#
    else :
            check1 = False


    #### check 2, is the object in the comments file and in the catalog?

    if (os.path.exists(catalog)) & (os.path.exists(comments)):

        ## read catalog data.
        #M.D.R# catalog_data =asciitable.read(catalog, format = 'sextractor')
        catalog_data =asciitable.read(catalog) #M.D.R#
        catalog_objid = catalog_data['ObjID']
        #print('\ncatalog_objid =', catalog_objid) #M.D.R#

        ### read comments file, get list of keeps vs. rejects
        f = open(comments)
        obj_comments = []
        reject_list = []

        for line in f:
           x = line.split()
           obj_comments.append(int(x[1]))

           if len(x) > 2 :
              if x[2] == 'rejected':
                    reject_list.append('rejected')
              else:
                    reject_list.append('keep')
           else:
                reject_list.append('keep')
        f.close()

    else:
         check2 = False
         check3 = False




    ####  is objid_parent in obj_comments and catalog_objid
    if (objid_parent[i] in  obj_comments) and (objid_parent[i] in catalog_objid):
        check2 = True
        #print('check2 passed!\n') #M.D.R#
        #print('objid_parent[i] =', objid_parent[i]) #M.D.R#
    else:
        check2 = False


    obj_comments = np.array(obj_comments)
    w=np.where(obj_comments == objid_parent[i])

    if np.size(w[0]) == 0 :
            check3 = False
    else:
        if reject_list[w[0][0]] == 'keep' :
            check3 = True
            #print('\ncheck3 passed!') #M.D.R#
        else:
            check3 = False

    #print par_parent[i], objid_parent[i], check1, check2, check3

    if (check1 == True) & (check2 == True) & (check3 == True) :
        keepers.append(1)
        #### GATHER INFORMATION FOR KEEPERS
        w=np.where(catalog_objid == objid_parent[i])
        z_wisp.append(catalog_data['redshift'][w[0][-1]])  ###  just in case the object is in the catlag more than once, take the last one.
        foiii_wisp.append(catalog_data['oiii_flux'][w[0][-1]])
        eoiii_wisp.append(catalog_data['oiii_err'][w[0][-1]])
        fha_wisp.append(catalog_data['hanii_flux'][w[0][-1]])
        eha_wisp.append(catalog_data['hanii_err'][w[0][-1]])
        fhb_wisp.append(catalog_data['hb_flux'][w[0][-1]])
        ehb_wisp.append(catalog_data['hb_err'][w[0][-1]])
        fsii_wisp.append(catalog_data['sii_flux'][w[0][-1]])
        esii_wisp.append(catalog_data['sii_err'][w[0][-1]])
        foii_wisp.append(catalog_data['oii_flux'][w[0][-1]])
        eoii_wisp.append(catalog_data['oii_err'][w[0][-1]])
        ew_oiii_obs_wisp.append(catalog_data['oiii_EW_obs'][w[0][-1]])

        print 'Keeper found! All checks passed for objid_parent[i]: ' + str(objid_parent[i]) #M.D.R#



    else:
        keepers.append(0)


z_wisp  = np.array(z_wisp)
foiii_wisp = np.array(foiii_wisp)
eoiii_wisp  = np.array(eoiii_wisp)
keepers = np.array(keepers)

print ''  #M.D.R#
print 'Size of Parent WISP catalog'
print np.size(masscat_wisp['ID'])
print ''  #M.D.R#
print 'Size of Retained WISP catalog'
print np.sum(keepers)



w=np.where(keepers == 1)

id_wisp = masscat_wisp['ID'][w]
#M.D.R# par_wisp = wispcat['par'][w]
par_wisp = 0 #M.D.R#
#M.D.R# logm_wisp = masscat_wisp['logM_50'][w]
logm_wisp = masscat_wisp['LMSTAR_50_SPPH'][w] #M.D.R#




#print np.size(z_wisp), np.size(id_wisp), np.size(par_wisp) #M.D.R#


#### let's make an output/master catalog of these things, so I can play around and decide on binning.

field = []

#M.D.R# for par in par_wisp:
#M.D.R#     field.append('Par' + str(par))
field.append('Par0') #M.D.R#
#M.D.R# for ff  in field_3dhst :
#M.D.R#     field.append(ff)

field = np.array(field)

#M.D.R# objid = np.append(id_wisp, id_3dhst)
#M.D.R# logm = np.append(logm_wisp, logm_3dhst)
#M.D.R# z = np.append(z_wisp, z_3dhst)
#M.D.R# foiii = np.append(foiii_wisp, foiii_3dhst)
#M.D.R# eoiii = np.append(eoiii_wisp, eoiii_3dhst)
#M.D.R# fha = np.append(fha_wisp, fha_3dhst)
#M.D.R# eha = np.append(eha_wisp, eha_3dhst)
#M.D.R# fhb = np.append(fhb_wisp, fhb_3dhst)
#M.D.R# ehb = np.append(ehb_wisp, ehb_3dhst)
#M.D.R# fsii = np.append(fsii_wisp, fsii_3dhst)
#M.D.R# esii = np.append(esii_wisp, esii_3dhst)
#M.D.R# foii = np.append(foii_wisp, foii_3dhst)
#M.D.R# eoii  = np.append(eoii_wisp, eoii_3dhst)
#M.D.R# ew_oiii_obs = np.append(ew_oiii_obs_wisp, ew_oiii_obs_3dhst)

# This was not done above for any lines except [O III]. #M.D.R#
# This is required in order for "w" to be indexed below.  #M.D.R#
fha_wisp = np.array(fha_wisp) #M.D.R#
eha_wisp  = np.array(eha_wisp) #M.D.R#
fhb_wisp = np.array(fhb_wisp) #M.D.R#
ehb_wisp  = np.array(ehb_wisp) #M.D.R#
fsii_wisp = np.array(fsii_wisp) #M.D.R#
esii_wisp  = np.array(esii_wisp) #M.D.R#
foii_wisp = np.array(foii_wisp) #M.D.R#
eoii_wisp  = np.array(eoii_wisp) #M.D.R#
ew_oiii_obs_wisp  = np.array(ew_oiii_obs_wisp) #M.D.R#

objid = id_wisp #M.D.R#
logm = logm_wisp #M.D.R#
z = z_wisp #M.D.R#
foiii = foiii_wisp #M.D.R#
eoiii = eoiii_wisp #M.D.R#
fha = fha_wisp #M.D.R#
eha = eha_wisp #M.D.R#
fhb = fhb_wisp #M.D.R#
ehb = ehb_wisp #M.D.R#
fsii = fsii_wisp #M.D.R#
esii = esii_wisp #M.D.R#
foii = foii_wisp #M.D.R#
eoii  = eoii_wisp #M.D.R#
ew_oiii_obs = ew_oiii_obs_wisp #M.D.R#

#print 'foiii = foiii_wisp' #M.D.R#
#print foiii #M.D.R#
#print type(foiii) #M.D.R#
#print 'fha = fha_wisp' #M.D.R#
#print fha #M.D.R#
#print type(fha) #M.D.R#
#print '' #M.D.R#



### remove a few more bad objects with zero error or foiii = -1
#M.D.R# w=np.where( (foiii > 0) & (eoiii > 0) & (z > 1.2))

# Minimum criteria.
w=np.where((foiii > 0) & (eoiii > 0) & (fhb > 0) & (ehb > 0)) #M.D.R#
# Where oiii is > median.
w=np.where((foiii > 0) & (eoiii > 0) & (foiii / eoiii > 3)) #M.D.R#
#print '5x median is: ' + str(np.median(foiii)*5.0) #M.D.R#
print '' #M.D.R#
#w=np.where((foiii > 0) & (eoiii > 0)) #M.D.R#
#w_tup=np.where( (foiii > 0) & (eoiii > 0) & (z > 1.2)) #M.D.R#
#w_lst = [] #M.D.R#
#w_lst.append(w_tup) #M.D.R#
#w = [x for xs in w_lst for x in xs] #M.D.R#

#print 'w =' #M.D.R#
#print w #M.D.R#
#print type(w) #M.D.R#
#w = np.asarray(w) #M.D.R#
#w = w.flatten() #M.D.R#
#print w #M.D.R#
#print type(w) #M.D.R#
#print 'w[0] = ' + str(w[0]) #M.D.R#
#print 'w[1] = ' + str(w[1]) #M.D.R#
#print 'w[1:5] = ' + str(w[1:5]) #M.D.R#
#print 'w[:] = ' + str(w[:]) #M.D.R#

#print 'fha = ' + str(fha)
#print 'type(fha) ', type(fha)
#print 'fha[0] = ' + str(fha[0]) #M.D.R#
#print 'fha[3] = ' + str(fha[3]) #M.D.R#
#print 'fha[w[1:5]] = ' + str(fha[w[1:5]]) #M.D.R#
#print 'fha[w[:]] = ' + str(fha[w[:]]) #M.D.R#
#print 'fha[w] = ' + str(fha[w]) #M.D.R#

#print 'fhb = ' + str(fhb)
#print 'fhb[w] = ' + str(fhb[w])

#M.D.R# field = field[w]
field = field[0] #M.D.R#
objid = objid[w]
logm = logm[w]
z = z[w]
foiii = foiii[w]
eoiii = eoiii[w]
fha = fha[w]
eha = eha[w]
fhb =fhb[w]
ehb = ehb[w]
fsii = fsii[w]
esii = esii[w]
foii = foii[w]
eoii = eoii[w]
ew_oiii_obs = ew_oiii_obs[w]


#M.D.R# print 'number of objects after OIII no cov,  OIII err = 0, z<1.2  removed'
print 'number of objects after clean sweep' #M.D.R# \
print np.size(field), np.size(objid), np.size(logm), np.size(z), np.size(foiii), np.size(eoiii)

objid = objid.astype(int)

#mex_agn_selection
ehb2 = copy.deepcopy(ehb)
w=np.where(fhb == 0 )
ehb2[w] = eoiii[w]


mex_flag = np.zeros(len(z))
hbsnr = fhb / ehb2
oiiihb = np.log10(foiii/fhb)
w=np.where(hbsnr <=3)
oiiihb[w] =  np.log10(foiii[w] / (3 * ehb2[w]))

x =logm
x2 = x-0.75
ymex_lower_obj = 0.375 / (x2 - 10.5) + 1.14  + 0.11
w=np.where(x2 > 9.6)
ymex_lower_obj[w] = 352.066 - 93.8249 * x2[w] + 8.32651 * x2[w]**2 - 0.246416 * x2[w]**3   + 0.11
agn_flag = np.where((oiiihb > ymex_lower_obj) & (logm >0))
mex_flag[agn_flag] = 1


xray_agn_flag = np.zeros(len(z))

#M.D.R# xray_field = xray['col1']
#M.D.R# xray_id = xray['col2']
xray_field = np.zeros(len(z)) #M.D.R# NOT USING.
xray_id = np.zeros(len(z)) #M.D.R# NOT USING.

for i in np.arange(len(xray_id)):
    w=np.where( (field == xray_field[i]) & (objid == xray_id[i]))
    if len(w[0]) == 1:
        xray_agn_flag[w] = 1.

ir_agn_flag = np.zeros(len(z))
#M.D.R# ir_field = donley['col1']
#M.D.R# ir_id = donley['col2']
ir_field = np.zeros(len(z)) #M.D.R#
ir_id = np.zeros(len(z)) #M.D.R#


for i in np.arange(len(ir_id)):
    w=np.where((field == ir_field[i]) & (objid == ir_id[i]))
    if len(w[0]) == 1:
        ir_agn_flag[w] = 1.



field = ['Par0'] * len(z) #M.D.R#
data = [field, objid, logm, z, foiii, eoiii, ew_oiii_obs, fhb, ehb, fha, eha, fsii, esii, foii, eoii, mex_flag, xray_agn_flag, ir_agn_flag]
#M.D.R# print data
colnames = ['Field', 'ID', 'logM', 'z', 'foiii', 'eoiii', 'EW_oiii_obs', 'fhb', 'ehb', 'fhanii', 'ehanii', 'fsii', 'esii', 'foii', 'eoii', 'MeX_AGN', 'Xray', 'IR_AGN']

asciitable.write(data, 'master_stacklist_run07_i02.dat', names = colnames, overwrite = True)
