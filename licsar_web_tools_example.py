#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 12:41:47 2021

@author: matthew
"""

import sys
import pdb

#%% Imports

from licsar_web_tools_downloading import download_LiCSAR_portal_data
from licsar_web_tools_converting import convert_LiCSAR_portal_data

import matplotlib.pyplot as plt

#%% An example on the North Anatolian fault

frameID = '014A_04939_131313'                                       # As defined on the LiCSAR portal
date_start = 20190901                                               # YYYYMMDD
date_end   = 20200101                                               # YYYYMMDD
download_metadata = True                                            # the network of baselines, DEM etc.  
between_epoch_files = ['geo.unw.png', 'geo.cc.png']                 # Possible files: 'geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif'
epoch_files = ['geo.mli.png','ztd.jpg']                             # Possible files: 'geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg'
n_para = 4                                                          # Parallelisation.  The number of cores is a good starting point.    


download_LiCSAR_portal_data(frameID, date_start, date_end, download_metadata, epoch_files, between_epoch_files, n_para)

sys.exit()

#%% An example of Sierra Negra volcano (Galapagos Archipelago, Ecuador)


frameID = '128D_09016_110500'                         #  As defined on the LiCSAR portal
date_start = 20180301                                 # YYYYMMDD
date_end   = 20180601                                 # YYYYMMDD
download_metadata = False                             # the network of baselines, DEM etc.  
between_epoch_files = ['geo.unw.png']                 # Possible files: 'geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif'
epoch_files = None                                    # Possible files: 'geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg'
n_para = 4                                            # Parallelisation.  The number of cores is a good starting point.    

download_LiCSAR_portal_data(frameID, date_start, date_end, download_metadata, epoch_files, between_epoch_files, n_para)

    
    
    
    
#%% Downloading some tif data

frameID = '162A_06192_060402'                         #  As defined on the LiCSAR portal
date_start = 20190101                                 # YYYYMMDD
date_end   = 20190201                                 # YYYYMMDD
download_metadata = True                              # the network of baselines, DEM etc.  
between_epoch_files = ['geo.unw.png', 'geo.unw.tif']  # Possible files: 'geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif'
epoch_files = ['geo.mli.tif']                         # Possible files: 'geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg'
n_para = 4                                            # Parallelisation.  The number of cores is a good starting point.    

download_LiCSAR_portal_data(frameID, date_start, date_end, download_metadata, epoch_files, between_epoch_files, n_para)


#%% And show how this can be converted to a simple to use numpy array.  

convert_metadata = True
metadata_nps, epoch_files_nps, between_epoch_files_nps = convert_LiCSAR_portal_data(frameID, convert_metadata, epoch_files, between_epoch_files, mask_vals = [0])


date1 = list(epoch_files_nps['geo.mli.tif'].keys())[0]                                                             # get the date of one of the images
f, ax = plt.subplots(1,1)
ax.imshow(epoch_files_nps['geo.mli.tif'][date1])                                                                    # and plot that date


date1 = list(between_epoch_files_nps['geo.unw.tif'].keys())[0]                                                     # get the date of one of the images
f, ax = plt.subplots(1,1)
ax.imshow(between_epoch_files_nps['geo.unw.tif'][date1])                                                            # and plot that date







