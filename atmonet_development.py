#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 12:41:47 2021

@author: matthew
"""
dependency_paths = {'syinterferopy_bin'  : '/home/matthew/university_work/01_blind_signal_separation_python/SyInterferoPy/lib/',              # Available from Github: https://github.com/matthew-gaddes/SyInterferoPy
                    'srtm_dem_tools_bin' : '/home/matthew/university_work/11_DEM_tools/SRTM-DEM-tools/',                                      # Available from Github: https://github.com/matthew-gaddes/SRTM-DEM-tools
                    'local_scripts'      : '/home/matthew/university_work/python_stuff/python_scripts/',
                    'LiCSBAS_1.3.9'      : '/home/matthew/university_work/03_automatic_detection_algorithm/07_LiCSBAS/LiCSBAS-1.3.9/bin/'}       #

import sys
for dependency_name, dependency_path in dependency_paths.items():
    sys.path.append(dependency_path)

#%% Imports

import rasterio
from pathlib import Path
import numpy as np
import glob

from small_plot_functions import matrix_show
# #%% Things to set
# frame = '074A_05055_131313'



# #%%

# ifgs = glob.glob(f"./{frame}/interferograms/*/*")

# for ifg in ifgs:
#     bands = rasterio.open(ifg)
#     unw_ph = bands.read()
    
    
    
    
# band_ds = rasterio.open(band_file)                                      # open the Geotiff
# if band_n == 0:                                                         # if it's the first one... 
#     ny, nx = band_ds.read(1).shape                                      # use it to get the size of the data
#     image_data = np.zeros((ny, nx, n_bands))                            # and initiate an array to store the data in,  channels last format.  Note - this won't work if channels with different resolution are used.  
# image_data[:,:, band_n] = band_ds.read(1)                               # read the geotif and add to the array that stores all the data 


#%% Gacos

# import os

# wd = os.getcwd()
# LiCSARweb = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'

# gacosdir = os.path.join(wd, 'GACOS')
# if not os.path.exists(gacosdir): os.mkdir(gacosdir)

# ### Get available dates
# print('\nDownload GACOS data', flush=True)
# url = os.path.join(LiCSARweb, trackID, frameID, 'epochs')
# response = requests.get(url)
# response.encoding = response.apparent_encoding #avoid garble
# html_doc = response.text
# soup = BeautifulSoup(html_doc, "html.parser")
# tags = soup.find_all(href=re.compile(r"\d{8}"))
# imdates_all = [tag.get("href")[0:8] for tag in tags]
# _imdates = np.int32(np.array(imdates_all))
# _imdates = (_imdates[(_imdates>=startdate)*(_imdates<=enddate)]).astype('str').tolist()
# print('  There are {} epochs from {} to {}'.format(len(_imdates),
#                                startdate, enddate), flush=True)

# ### Extract available dates
# print('  Searching available epochs ({} parallel)...'.format(n_para), flush=True)

# args = [(i, len(_imdates),
#          os.path.join(url, imd, '{}.sltd.geo.tif'.format(imd)),
#          os.path.join(gacosdir, imd+'.sltd.geo.tif')
#          ) for i, imd in enumerate(_imdates)]

# p = multi.Pool(n_para)
# rc = p.map(check_gacos_wrapper, args)
# p.close()

# n_im_existing = 0
# n_im_unavailable = 0
# imdates_dl = []
# for i, rc1 in enumerate(rc):
#     if rc1 == 0:  ## No need to download
#         n_im_existing = n_im_existing + 1
#     if rc1 == 3 or rc1 == 5:  ## Can not download
#         n_im_unavailable = n_im_unavailable + 1
#     elif rc1 == 1 or rc1 == 2  or rc1 == 4:  ## Need download
#         imdates_dl.append(_imdates[i])

# n_im_dl = len(imdates_dl)

# if n_im_existing > 0:
#     print('  {} GACOS data already downloaded'.format(n_im_existing), flush=True)
# if n_im_unavailable > 0:
#     print('  {} GACOS data unavailable'.format(n_im_unavailable), flush=True)

# ### Download
# if n_im_dl > 0:
#     print('{} GACOS data will be downloaded'.format(n_im_dl), flush=True)
#     print('Download GACOS ({} parallel)...'.format(n_para), flush=True)
#     ### Download
#     args = [(i, imd, n_im_dl,
#              os.path.join(url, imd, '{}.sltd.geo.tif'.format(imd)),
#              os.path.join(gacosdir, '{}.sltd.geo.tif'.format(imd))
#              ) for i, imd in enumerate(imdates_dl)]
    
#     p = multi.Pool(n_para)
#     p.map(download_wrapper, args)
#     p.close()
# else:
#     print('No GACOS data available from {} to {}'.format(startdate, enddate), flush=True)


#%%

# import subprocess
# frameID = '074A_05055_131313'
# date_start = 20190901
# date_end = 20200101

# # #LiCSBAS01_get_geotiff.py [-f frameID] [-s yyyymmdd] [-e yyyymmdd] [--get_gacos] [--n_para int]

# subprocess.call(f"LiCSBAS01_get_geotiff.py -f {frameID} -s {date_start} -e {date_end} --get_gacos {True}" , shell=True)                 # 

# # #subprocess.call(f"LiCSBAS02_ml_prep.py -i {GEOCdir} -o {GEOCmldir} -n {downsampling} --n_para {n_para}" + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)                 # This creates the files in GEOCmlXXX, including the png preview of unw, note that 1 is stdout, -a to append

# import sys; sys.exit()

#%%


#ver=1.6; date=20200911; author="Y. Morishita"
# print("\n{} ver{} {} {}".format(os.path.basename(argv[0]), ver, date, author), flush=True)
# print("{} {}".format(os.path.basename(argv[0]), ' '.join(argv[1:])), flush=True)

#%% Dependencies

from pathlib import Path

import os
import datetime as dt
import requests
from bs4 import BeautifulSoup                                                                            # to get data out of HTML etc.  
import re
import multiprocessing as multi

def comp_size_time(file_remote, file_local):
    """
    Compare size and time of remote and local files.
    Returns:
        0 : Size is identical and local time is new
        1 : Size is not identical
        2 : Size is identical but remote time is new
        3 : Remote not exist
    """
    import requests
    import datetime as dt
    import dateutil

    response = requests.head(file_remote, allow_redirects=True)
    
    if response.status_code != 200:
        return 3

    size_remote = int(response.headers.get("Content-Length"))
    size_local = os.path.getsize(file_local)
    if size_remote != size_local: ### Different size
        return 1

    time_remote = dateutil.parser.parse(response.headers.get("Last-Modified"))
    time_local = dt.datetime.fromtimestamp(os.path.getmtime(file_local), dt.timezone.utc)
    if time_remote > time_local: ### New file
        return 2

    return 0

def download_data(url, file):
    """ Download a file from a url and set the name of the downloaded file.  
    Inputs:
        url | Path | url of file to download (i.e. including the file name)
        file | string | name of file to save as
    Returns:
        file
    History:
        2021_??_?? | YM | Part of LiCSBAS
        2021_03_25 | MEG | Write the docs.   
        
    """
    import time
    import requests
    
    def convert_size(size_bytes):
        if size_bytes == 0:
            return "0B"
        size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
        i = int(np.floor(np.log(size_bytes)/np.log(1024)))
        p = np.power(1024, i)
        s = round(size_bytes / p, 2)
        return "%s%s" % (s, size_name[i])

    try:
       start = time.time() 
       with requests.get(url) as res:
            res.raise_for_status()
            with open(file, "wb") as output:
                output.write(res.content)
       elapsed = int(time.time()-start)
       if elapsed==0: elapsed=elapsed+1
       fsize = convert_size(os.path.getsize(file))
       speed = convert_size(int(os.path.getsize(file)/elapsed))
       print('  {}, {}, {}s, {}/s'.format(os.path.basename(file),
                                          fsize, elapsed, speed))

    except:
        print(
            "    Error while downloading from {}".format(url),
            file=sys.stderr,
            flush=True,
        )
        return


def download_wrapper(args):
    i, ifgd, n_dl, url_data, path_data = args
    dir_data = os.path.dirname(path_data)
    print('  Donwnloading {} ({}/{})...'.format(ifgd, i+1, n_dl), flush=True)
    if not os.path.exists(dir_data): os.mkdir(dir_data)
    download_data(url_data, path_data)

def check_gacos_wrapper(args):
    """
    Returns :
        0 : Local exist, complete, and new (no need to donwload)
        1 : Local incomplete (need to re-donwload)
        2 : Local old (no need to re-donwload)
        3 : Remote not exist  (can not compare, no download)
        4 : Local not exist and remote exist (need to download)
        5 : Local not exist but remote not exist (can not download)
    """
    i, n_data, url_data, path_data = args
    bname_data = os.path.basename(path_data)
    
    if np.mod(i, 10) == 0:
        print("  {0:3}/{1:3}".format(i, n_data), flush=True)

    if os.path.exists(path_data):
        rc = comp_size_time(url_data, path_data)
        if rc == 1:
            print("Size of {} is not identical.".format(bname_data), flush=True)
        elif rc == 2:
            print("Newer {} available.".format(bname_data), flush=True)
        return rc
    else:
        response = requests.head(url_data, allow_redirects=True)
        if response.status_code == 200:
            return 4
        else:
            return 5
   
    

def ifgdates2imdates(ifgdates):
    primarylist = []
    secondarylist = []
    for ifgd in ifgdates:
        primarylist.append(ifgd[:8])
        secondarylist.append(ifgd[-8:])

    imdates = list(set(primarylist+secondarylist)) # set is a unique operator
    imdates.sort()

    return imdates

def check_exist_wrapper(args):
    """
    Returns :
        0 : Local exist, complete, and new (no need to donwload)
        1 : Local incomplete (need to re-donwload)
        2 : Local old (no need to re-donwload)
        3 : Remote not exist  (can not compare, no download)
        4 : Local not exist (need to download)
    """

    i, n_data, url_data, path_data = args
    bname_data = os.path.basename(path_data)
    
#    if np.mod(i, 10) == 0:
#        print("  {0:3}/{1:3}".format(i, n_data), flush=True)

    if os.path.exists(path_data):
        rc = comp_size_time(url_data, path_data)
        if rc == 1:
            print("Size of {} is not identical.".format(bname_data), flush=True)
        elif rc == 2:
            print("Newer {} available.".format(bname_data), flush=True)
        return rc
    else:
        return 4
    


#%% Set default
frameID = '074A_05055_131313'                                                                   # Central Spain
frameID = '014A_04939_131313'                                                                   # North Anatolian Fault
startdate = 20141001
enddate = int(dt.date.today().strftime("%Y%m%d"))
get_gacos = True
n_para = 4

# it expects you to run the file in te frames dir.  I ahven't done that so it's all a mess.  

#%% Determine frameID
# wd = os.getcwd()
# #wd = f"./{frameID}/"
# if not os.path.exists(wd): os.mkdir(wd)

trackID = str(int(frameID[0:3]))


#%% Directory and file setting
# outdir = os.path.join(wd, 'GEOC')
# if not os.path.exists(outdir): os.mkdir(outdir)
# os.chdir(outdir)

LiCSARweb = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'


#%% Download the ENU and hgt files

def download_aux_file(LiCSARweb, trackID, frameId):
    """
    """
    #import pdb; pdb.set_trace()

    for ENU in ['E', 'N', 'U', 'hgt']:
        enutif = f"{frameID}.geo.{ENU}.tif"                                                                 # build a string of the name of the file
        url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', enutif)                                 # build the url to where that file would be stored by LiCSARweb
        if os.path.exists(enutif):                                                                          # see if the file exists locally
            rc = comp_size_time(url, enutif)                                                                # if it does, compare the remote (first) and local (second) file.  
            if rc == 0:
                print('{} already exist. Skip download.'.format(enutif), flush=True)
                continue                                                                                    # jumps to next part of for loop so skips download
            elif rc == 3:
                print('{} not available. Skip download.'.format(enutif), flush=True)
                continue                                                                                    # jumps to next part of for loop so skips download
            else:
                if rc == 1:
                    print("Size of {} is not identical.".format(enutif))                                    # download will occur
                elif rc == 2:
                    print("Newer {} available.".format(enutif))                                             # download will occur
        
        print('Download {}'.format(enutif), flush=True)                                                     # update user
        download_data(url, enutif)                                                                          # and download.  enturif is path for download.  
    
    
    
    print('Download baselines', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'baselines')                                # construct the URL for the baselines file
    download_data(url, 'baselines')                                                                         # download the baselines file.  Of the form YYYYMMDD YYYYMMDD perp_baseline temporal_baseline
    
    print('Download metadata.txt', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'metadata.txt')                             # construct the URL for the metadata file
    download_data(url, 'metadata.txt')                                                                      # master date and frame centre acquisition time.  
    
    print('Download network.png', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'network.png')                             # construct the URL for the metadata file
    download_data(url, 'network.png')                                                                      # master date and frame centre acquisition time.  

download_aux_file(LiCSARweb, trackID, frameID)

import sys; sys.exit()

#%% mli
# mlitif = frameID+'.geo.mli.tif'
# if os.path.exists(mlitif):
#     print('{} already exist. Skip.'.format(mlitif), flush=True)
# else:
#     ### Get available dates
#     print('Searching earliest epoch for mli...', flush=True)
#     url = os.path.join(LiCSARweb, trackID, frameID, 'epochs')
#     response = requests.get(url)
    
#     response.encoding = response.apparent_encoding #avoid garble
#     html_doc = response.text
#     soup = BeautifulSoup(html_doc, "html.parser")
#     tags = soup.find_all(href=re.compile(r"\d{8}"))
#     imdates_all = [tag.get("href")[0:8] for tag in tags]
#     _imdates = np.int32(np.array(imdates_all))
#     _imdates = (_imdates[(_imdates>=startdate)*(_imdates<=enddate)]).astype('str').tolist()
    
#     ## Find earliest date in which mli is available
#     imd1 = []
#     for i, imd in enumerate(_imdates):
#         if np.mod(i, 10) == 0:
#             print("\r  {0:3}/{1:3}".format(i, len(_imdates)), end='', flush=True)
#         url_epoch = os.path.join(url, imd)
#         response = requests.get(url_epoch)
#         response.encoding = response.apparent_encoding #avoid garble
#         html_doc = response.text
#         soup = BeautifulSoup(html_doc, "html.parser")
#         tag = soup.find(href=re.compile(r"\d{8}.geo.mli.tif"))
#         if tag is not None:
#             print('\n{} found as earliest.'.format(imd))
#             imd1 = imd
#             break

#     ### Download
#     if imd1:
#         print('Donwnloading {}.geo.mli.tif as {}.geo.mli.tif...'.format(imd1, frameID), flush=True)
#         url_mli = os.path.join(url, imd1, imd1+'.geo.mli.tif')
#         tools_lib.download_data(url_mli, mlitif)
#     else:
#         print('\nNo mli available on {}'.format(url), file=sys.stderr, flush=True)


#%% GACOS if specified
os.chdir('..')
if get_gacos:
    gacosdir = os.path.join(wd, 'GACOS')                                                                        # name of directory for GACOS products
    if not os.path.exists(gacosdir): os.mkdir(gacosdir)                                                         # if it doesn't exist, make it.  

    ### Get available dates
    print('\nDownload GACOS data', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'epochs')
    response = requests.get(url)                                                                                # get the webpage with all the epochs on
    response.encoding = response.apparent_encoding                                                              # avoid garble?
    html_doc = response.text                                                                                    # turns the text into a giant string.  
    soup = BeautifulSoup(html_doc, "html.parser")                                                               # break the string into something that looks more like normal html
    tags = soup.find_all(href=re.compile(r"\d{8}"))                                                             # get just the links (which are epochs) from the html
    imdates_all = [tag.get("href")[0:8] for tag in tags]                                                        # get the epochs (as a list)
    _imdates = np.int32(np.array(imdates_all))                                                                  # convert to numpy array
    _imdates = (_imdates[(_imdates>=startdate)*(_imdates<=enddate)]).astype('str').tolist()                     # back to list?
    print('  There are {} epochs from {} to {}'.format(len(_imdates),                                           # len(_imdates) is the number of epochs.  
                                   startdate, enddate), flush=True)

    ### Extract available dates
    print('  Searching available epochs ({} parallel)...'.format(n_para), flush=True)

    args = [(i, len(_imdates),                                                                                  # a list with an entry for each epoch.  
             os.path.join(url, imd, '{}.sltd.geo.tif'.format(imd)),                                             # each item in the list is a tuple, with 4 entries.  
             os.path.join(gacosdir, imd+'.sltd.geo.tif')                                                         # the first is?  The second is?  
             ) for i, imd in enumerate(_imdates)]                                                               # the 3rd is the path to the sltd file from Gacod, and the 4th is the path so save that to.  

    p = multi.Pool(n_para)
    rc = p.map(check_gacos_wrapper, args)                                                                       # for all epochs, check if Gacos data needs to be downloaded (compare locate [if exists] and remote)
    p.close()                                                                                                   # rc is as list of the same length as the number of epochs.  

    n_im_existing = 0
    n_im_unavailable = 0
    imdates_dl = []
    for i, rc1 in enumerate(rc):                                                                                 # for each epoch's rc value, determine if we will download or not
        if rc1 == 0:                                                                                             # No need to download
            n_im_existing = n_im_existing + 1
        if rc1 == 3 or rc1 == 5:                                                                                 # Cannot download
            n_im_unavailable = n_im_unavailable + 1
        elif rc1 == 1 or rc1 == 2  or rc1 == 4:                                                                  # Need to download
            imdates_dl.append(_imdates[i])                                                                       # append to a list of which dates to download

    n_im_dl = len(imdates_dl)

    if n_im_existing > 0:
        print('  {} GACOS data already downloaded'.format(n_im_existing), flush=True)
    if n_im_unavailable > 0:
        print('  {} GACOS data unavailable'.format(n_im_unavailable), flush=True)

    ### Download
    if n_im_dl > 0:
        print('{} GACOS data will be downloaded'.format(n_im_dl), flush=True)
        print('Download GACOS ({} parallel)...'.format(n_para), flush=True)
        ### Download
        args = [(i, imd, n_im_dl,
                 os.path.join(url, imd, '{}.sltd.geo.tif'.format(imd)),
                 os.path.join(gacosdir, '{}.sltd.geo.tif'.format(imd))
                 ) for i, imd in enumerate(imdates_dl)]
        
        p = multi.Pool(n_para)
        p.map(download_wrapper, args)
        p.close()
    else:
        print('No GACOS data available from {} to {}'.format(startdate, enddate), flush=True)


#%% unw and cc
### Get available dates
print('\nDownload geotiff of unw and cc', flush=True)
url_ifgdir = os.path.join(LiCSARweb, trackID, frameID, 'interferograms')
response = requests.get(url_ifgdir)

response.encoding = response.apparent_encoding #avoid garble
html_doc = response.text
soup = BeautifulSoup(html_doc, "html.parser")
tags = soup.find_all(href=re.compile(r"\d{8}_\d{8}"))
ifgdates_all = [tag.get("href")[0:17] for tag in tags]

### Extract during start_date to end_date
ifgdates = []
for ifgd in ifgdates_all:
    mimd = int(ifgd[:8])
    simd = int(ifgd[-8:])
    if mimd >= startdate and simd <= enddate:
        ifgdates.append(ifgd)

n_ifg = len(ifgdates)
imdates = ifgdates2imdates(ifgdates)
print('{} IFGs available from {} to {}'.format(n_ifg, imdates[0], imdates[-1]), flush=True)

### Check if both unw and cc already donwloaded, new, and same size
print('Checking existing unw and cc ({} parallel, may take time)...'.format(n_para), flush=True)

## unw
args = [(i, n_ifg,
         os.path.join(url_ifgdir, ifgd, '{}.geo.unw.tif'.format(ifgd)),
         os.path.join(ifgd, '{}.geo.unw.tif'.format(ifgd))
         ) for i, ifgd in enumerate(ifgdates)]

p = multi.Pool(n_para)
rc = p.map(check_exist_wrapper, args)
p.close()

n_unw_existing = 0
unwdates_dl = []
for i, rc1 in enumerate(rc):
    if rc1 == 0:  ## No need to download
        n_unw_existing = n_unw_existing + 1
    if rc1 == 3 or rc1 == 5:  ## Can not download
        print('  {}.geo.unw.tif not available.'.format(ifgdates[i]), flush=True)
    elif rc1 == 1 or rc1 == 2  or rc1 == 4:  ## Need download
        unwdates_dl.append(ifgdates[i])

## cc
args = [(i, n_ifg,
         os.path.join(url_ifgdir, ifgd, '{}.geo.cc.tif'.format(ifgd)),
         os.path.join(ifgd, '{}.geo.cc.tif'.format(ifgd))
         ) for i, ifgd in enumerate(ifgdates)]

p = multi.Pool(n_para)
rc = p.map(check_exist_wrapper, args)
p.close()

n_cc_existing = 0
ccdates_dl = []
for i, rc1 in enumerate(rc):
    if rc1 == 0:  ## No need to download
        n_cc_existing = n_cc_existing + 1
    if rc1 == 3 or rc1 == 5:  ## Can not download
        print('  {}.geo.cc.tif not available.'.format(ifgdates[i]), flush=True)
    elif rc1 == 1 or rc1 == 2  or rc1 == 4:  ## Need download
        ccdates_dl.append(ifgdates[i])

n_unw_dl = len(unwdates_dl)
n_cc_dl = len(ccdates_dl)
print('{} unw already downloaded'.format(n_unw_existing), flush=True)
print('{} unw will be downloaded'.format(n_unw_dl), flush=True)
print('{} cc already downloaded'.format(n_cc_existing), flush=True)
print('{} cc will be downloaded'.format(n_cc_dl), flush=True)

### Download unw with parallel
if n_unw_dl != 0:
    print('Download unw ({} parallel)...'.format(n_para), flush=True)
    args = [(i, ifgd, n_unw_dl,
             os.path.join(url_ifgdir, ifgd, '{}.geo.unw.tif'.format(ifgd)),
             os.path.join(ifgd, '{}.geo.unw.tif'.format(ifgd))
             ) for i, ifgd in enumerate(unwdates_dl)]
    
    p = multi.Pool(n_para)
    p.map(download_wrapper, args)
    p.close()
   
### Download cc with parallel
if n_cc_dl != 0:
    print('Download cc ({} parallel)...'.format(n_para), flush=True)
    args = [(i, ifgd, n_cc_dl,
             os.path.join(url_ifgdir, ifgd, '{}.geo.cc.tif'.format(ifgd)),
             os.path.join(ifgd, '{}.geo.cc.tif'.format(ifgd))
             ) for i, ifgd in enumerate(ccdates_dl)]
    
    p = multi.Pool(n_para)
    p.map(download_wrapper, args)
    p.close()
   

#%% Finish
# elapsed_time = time.time()-start
# hour = int(elapsed_time/3600)
# minite = int(np.mod((elapsed_time/60),60))
# sec = int(np.mod(elapsed_time,60))
# print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))

# print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
# print('Output directory: {}\n'.format(outdir))

