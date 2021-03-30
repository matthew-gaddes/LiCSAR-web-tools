#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 17:26:27 2021

@author: matthew
"""

#%%


def convert_LiCSAR_portal_data(frameID, convert_metadata, epoch_files, between_epoch_files, mask_vals = [0]):
    """Given a set of directories downloaded from the COMET-LiCS portal, convert them to numpy arrays.  
    Inputs:
        frameID              | string | The frameID, as defined on the COMET portal.  
        convert_metadata    | boolean | If True, the files in 'metadata' for that frame are converted.  E.g. the network of baselines, or DEM.  
        epoch_files          | None or list | List of strings with the extensions of the files that we wish to download.  There are per epoch (i.e things like radar intensity)  
                                              E.g. ['geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg']
        between_epoch_files | None or list | List of strings with the extensions of the files that we wish to download.  There are between epochs (e.g. phase)   
                                             ['geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif']:
        mask_vals | list | 0 is usually used for water/incoherent pixels, but this (or multiple values) can be set here.  

    Returns:
        metadata_nps            | dict |  E, N, U and hgt (dem) files each as numpy arrays.  
        epoch_files_nps          | dict of dicts | each type of file extension is a dictionary, and each one each date is one of the filews as a numpy array.  
        between_epoch_files_nps  | dict of dicts | each type of file extension is a dictionary, and each one each date is one of the filews as a numpy array.  
        
    History:
        2021_03_30 | MEG | Written
    
        
    """
    from pathlib import Path
    import glob
    import os

    
    # 0: Possibly conver the metadata tif files
    if convert_metadata:
        metadata_tifs = glob.glob(os.path.join(frameID, 'metadata', '*.tif'))                                       # get the paths to the three metadata tifs (E,N, U and hgt [the DEM])
        metadata_nps = licsar_tif_to_numpy(metadata_tifs, name_block = 2, mask_vals = mask_vals)                    # conver them to numpy arrays.                      
        
    
    # 1: Possibly convert the files that are processed per epoch
    epoch_files_nps = {}                                                                                                        # initiate to store each output.  
    if epoch_files is not None:
        for epoch_file in epoch_files:
            if epoch_file.split('.')[-1] == 'tif':                                                                              # check that the file to be converted is a tif
                epoch_file_tifs = glob.glob(os.path.join(frameID, 'epochs', '*', f'*{epoch_file}'))                             # get the path to all the files with that extension
                epoch_files_nps[epoch_file] = licsar_tif_to_numpy(epoch_file_tifs, name_block = 0, mask_vals = mask_vals)       # convert to a numpy array                                
            else:
                print(f"Files with the extension {epoch_file} aren't tifs so will be skipped.  ")                               # advise user if skipping as not a tif

    # 2: Possibly convert the files that are processed between epochs
    between_epoch_files_nps = {}                                                                                                                        # initiate to store each output.  
    if between_epoch_files is not None:
        for between_epoch_file in between_epoch_files:                                                                                                  # loop through each file extension
            if between_epoch_file.split('.')[-1] == 'tif':                                                                                              # check that it's a tif file
                between_epoch_file_tifs = glob.glob(os.path.join(frameID, 'interferograms', '*', f'*{between_epoch_file}'))                             # get a list of all the files.  
                between_epoch_files_nps[between_epoch_file] = licsar_tif_to_numpy(between_epoch_file_tifs, name_block = 0, mask_vals = mask_vals)       # conver them to a numpy array.                                          
            else:
                print(f"Files with the extension {between_epoch_file} aren't tifs so will be skipped.  ")                                               # advice user if skipping as not a tif.  
        
    return metadata_nps, epoch_files_nps, between_epoch_files_nps



#%%

def licsar_tif_to_numpy(licsar_tifs, name_block = 0, mask_vals = [0]):
    """ Gievn a list of paths to licsar tifs, return them as a dictionary of numpy array.  
    
    Inputs:
        licsar_tifs | list of paths (or strings) | paths to each tif file that is to be converted.  
        name_block | which part of the filename is used to name the numpy array, after the name is split by . (dot).  e.g. -1 would select .tff if it's geo.unw.tif
        mask_vals | list | 0 is usually used for water/incoherent pixels, but this (or multiple values) can be set here.  
    
    Returns:
        licsar_nps | dict | tiff data in numpy form.  
        
    History:
        2021_03_30 | MEG | Written
    """
    import numpy as np
    import gdal
    from pathlib import Path
    
    licsar_nps = {}                                                                                         # initiate a dictionary to store the numpy arrays in.  
    for licsar_tif in licsar_tifs:                                                                          # loop through each of the files.  
        licsar_np = np.float32(gdal.Open(licsar_tif).ReadAsArray())                                         # convert the tif to a numpy array.  
        for mask_val in mask_vals:
            licsar_np[licsar_np == mask_val] = np.nan                                                       # mask areas that are exactly equal to the masking value (which we're looping through)
        np_name = Path(licsar_tif).parts[-1].split('.')[name_block]                                         # break the path up to get just the filename, then break that on dots and select which part (set by name_block)
        licsar_nps[np_name] = licsar_np
            
    return licsar_nps
