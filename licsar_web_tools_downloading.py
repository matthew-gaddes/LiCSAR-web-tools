#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 11:13:12 2021

@author: matthew
"""

#%%

def download_LiCSAR_portal_data(frameID, date_start, date_end, download_metadata = True, epoch_files = None, between_epoch_files = None,
                                n_para = 1):
    """A function to download data from the COMET LiCS Sentinel-1 portal.  Uses many functions and ideas from Step-01 of LiCSBAS.  
    Inputs:
        frameID              | string | The frameID, as defined on the COMET portal.  
        date_start           | string | YYYYMMDD
        date_end             | string | YYYYMMDD
        download_metadata    | boolean | If True, the files in 'metadata' for that frame are downloaded.  E.g. the network of baselines, or DEM.  
        epoch_files          | None or list | List of strings with the extensions of the files that we wish to download.  There are per epoch (i.e things like radar intensity)  
                                              E.g. ['geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg']
        between_epoch_files | None or list | List of strings with the extensions of the files that we wish to download.  There are between epochs (e.g. phase)   
                                             ['geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif']:
         n_para             | int | Sets the parallelisation of downloads.  Number of CPU cores is a good starting point.  
     Returns:
         Local files.  
     History:
         2021_03_29 | MEG | Written.  
                                                    
    """
    
    
    LiCSARweb = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'                                      # Shouldn't need to change
    
    import os
    from pathlib import Path
    
    trackID = str(int(frameID[0:3]))
    
    print(f"\n\n\n\nStarting LiCSAR_web_tools for frame {frameID}.")
    # -1: Check inputs for common mistakes (not exhaustive)
    if between_epoch_files is not None:
        for between_epoch_file in between_epoch_files:
            if between_epoch_file not in ['geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif']:
                print(f"'{between_epoch_file}' is not a LiCSAR file extension.  No attempt will be made to download it, but we'll continue anyway.  ")
                if len(between_epoch_files) == 1:
                    between_epoch_files = None
                else:
                    between_epoch_files.remove(between_epoch_file)
        
    if epoch_files is not None:
        for epoch_file in epoch_files:
            if epoch_file not in ['geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg']:
                print(f"'{epoch_file}' is not a LiCSAR file extension.  No attempt will be made to download it, but we'll continue anyway.  ")
                if len(epoch_files) == 1:
                    epoch_files = None
                else:
                    epoch_files.remove(epoch_file)
    print(f"Downloading of files with extensions {between_epoch_files} (between epochs) and {epoch_files} (per epoch) between {date_start} and {date_end} will be attempted.  ")            
    
    if int(date_start) > int(date_end):
        raise Exception(f"The start date ({date_start}) is after the end date ({date_end}).  Exiting.  ")
    
    
    
    # 0: Build the local directories for things to be saved in.  Mirror of the system used on LiCSAR hub
    directories_dict = {'frame'           : Path(f"./{frameID}")}                                              # this is the same struture as the COMET LiCS portal.  
    directories_dict = {'interferograms'  : directories_dict['frame'] / 'interferograms',
                        'epochs'          : directories_dict['frame'] / 'epochs',
                        'metadata'        : directories_dict['frame'] / 'metadata'}
                        
    for directory in directories_dict.values():
        if not os.path.exists(directory): os.makedirs(directory)                                            # note use of makedirs to required intermediate directoryies can be made if needed.  
        
    # 1: Possibly download the metadata files
    if download_metadata:
        download_LiCSAR_metadata_files(trackID, frameID, directories_dict['metadata'])
        
    # 2: Possibly download files that are processed per epoch (e.g. Gacos)
    if epoch_files is not None:
        url_epochdir = os.path.join(LiCSARweb, trackID, frameID, 'epochs')
        for epoch_file in epoch_files:
            download_LiCSAR_files(url_epochdir, directories_dict['epochs'], date_start, date_end, epoch_file, double_date=False, n_para=n_para)
    
    
    # 3: Possibly download files that are processed between epochs (e.g. interferograms)
    if between_epoch_files is not None:
        url_ifgdir = os.path.join(LiCSARweb, trackID, frameID, 'interferograms')
        for between_epoch_file in between_epoch_files:
            download_LiCSAR_files(url_ifgdir, directories_dict['interferograms'], date_start, date_end, between_epoch_file, n_para=n_para)
    



#%%


def download_LiCSAR_files(remote_dir, local_dir, date_start, date_end, file_extension, double_date = True, n_para = 1):
    """ A function to download  any type of file that is either processed per epoch (e.g. Gacos files), or between epochs (e.g. unwrapped phase.  )
    Modified from scripts contained in LiCSBAS.  
    
    Inputs:
        remote_dir | string | the remote directory with the dates for each folder in.  e.g. http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/14/014A_04939_131313/interferograms'
                              Note that this is handled as a string as it's not dependent on the users' OS (i.e. it's aways coming from a linux server)
        local_dir | string (or Path) |  the local directory with the dates for each folder in'014A_04939_131313/interferograms'.  Caution if using strings for paths!
        date_start | string | YYYYMMDD form
        date_end | string | YYYYMMDD form
        file_extension | string | e.g. 'geo.unw.tif'  No . at the start!
        double_date | boolean | set to True if working with interfeograms, but false if Epochs.  
        n_para | int | number of parallel processes possible. Number of cores is a good starting point.  
    Returns:
        Local files.  
    History:
        2021_03_26 | MEG | Created from modified LiCSBAS scripts.  
        
    """
    
    import os
    import multiprocessing as multi
    
    date_start = int(date_start)                                                                        # convert from a string to an integer
    date_end = int(date_end)
        
    # 0: Get the dates required:
    file_dates_all = get_folders_from_LiCSAR_portal(remote_dir, double_date = double_date)              # get all the dates available on the portal
    file_dates = []                                                                                     # initate for the dates we want
    for file_date in file_dates_all:
        if double_date:
            mimd = int(file_date[:8])
            simd = int(file_date[-8:])
            if mimd >= date_start and simd <= date_end:
                file_dates.append(file_date)
        else:
            imd = int(file_date[:8])
            if imd >= date_start and imd <= date_end:              
                file_dates.append(file_date)


    # 1:  Determine which files need to be downloaded
    n_files= len(file_dates)
    args = [(i, n_files,
             os.path.join(remote_dir, file_date, f"{file_date}.{file_extension}"),
             os.path.join(local_dir, file_date, f"{file_date}.{file_extension}")
             ) for i, file_date in enumerate(file_dates)]
    p = multi.Pool(n_para)
    rc = p.map(check_exist_wrapper, args)                                                               # compare local (if it exists) and remote.  
    p.close()

    n_files_existing = 0
    files_dl = []
    for i, rc1 in enumerate(rc):                                                                    # rc is a list of the status for each file.  Loop through it.  
        if rc1 == 0:                                                                                  ## No need to download
            n_files_existing += 1
        if rc1 == 3:                                                                                 ## Can not download
            print(f' {file_dates[i]}.{file_extension} not available.', flush=True)
        elif rc1 == 1 or rc1 == 2  or rc1 == 4:                                                      ## Need download
            files_dl.append(file_dates[i])                                                              # if it does need downloading, append to list


    n_files_dl = len(files_dl)                                                                          # total number of files that will need to be downloaded
    print(f'{n_files_existing} {file_extension} files already downloaded', flush=True)
    print(f'{n_files_dl} {file_extension} files will be downloaded', flush=True)
    
    # 2: Do the downloading:
    if n_files_dl != 0:
        print('Downloading files ({} parallel)...'.format(n_para), flush=True)
        args = [(i, file_dl, n_files_dl,                                                            # creates a list of tuples, each item in the list is one file to be downloaded.  tuple is file n, date, remote path (to download from) and local path (to download to)
                 os.path.join(remote_dir, file_dl, f"{file_dl}.{file_extension}"),
                 os.path.join(local_dir, file_dl, f"{file_dl}.{file_extension}")
                 ) for i, file_dl in enumerate(files_dl)]
        
        p = multi.Pool(n_para)                                                                          # initiate the parallel downloading
        p.map(download_wrapper, args)                                                                   # do the downloading
        p.close()                                                                                       # close the parallelisation
  
    
def download_wrapper(args):
    """A function to download files.  Used inside a funtion to parallelise.   """
    import os
    i, ifgd, n_dl, url_data, path_data = args
    dir_data = os.path.dirname(path_data)
    print('  Donwnloading {} ({}/{})...'.format(ifgd, i+1, n_dl), flush=True)
    if not os.path.exists(dir_data): os.mkdir(dir_data)
    download_data(url_data, path_data)
    
    
#%%

def download_LiCSAR_metadata_files(trackID, frameID, outdir, LiCSARweb = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'):
    """The LiCSAR portal contains a folder of metadata, including products such as the DEM, or a figure showing the baselines.  
    This function downloads them.  
    Inputs:
        trackID | int | the track ID.  Usually the first part of the frameID
        frameID | str | the frame ID, as used on the COMET LiCS portal.  
        outdir | str of path | folder to save the outputs to
        LiCSARweb | string | parent directory for the COMET LiCS data.  Careful with how handled as it looks like a path, but it's really a URL.  
    Returns:
        E, N, U files.  
        hgt file                DEM
        baselines
        metadata
        network.png
    History:
        2021_03_26 | MEG | Modified from a LiCSBAS script.  
    """
    import os  

    for ENU in ['E', 'N', 'U', 'hgt']:
        enutif = f"{frameID}.geo.{ENU}.tif"                                                                 # build a string of the name of the file
        url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', enutif)                                 # build the url to where that file would be stored by LiCSARweb
        enutif = outdir / enutif                                                                            # conert from filename to path to file.  
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
        download_data(url, enutif)                                                                          # and download.  enturif is path for download, already converted to a posipath
    
    print('Download baselines', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'baselines')                                # construct the URL for the baselines file
    download_data(url, outdir / 'baselines')                                                                # download the baselines file.  Of the form YYYYMMDD YYYYMMDD perp_baseline temporal_baseline
    
    print('Download metadata.txt', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'metadata.txt')                             # construct the URL for the metadata file
    download_data(url, outdir /  'metadata.txt')                                                            # master date and frame centre acquisition time.  
    
    print('Download network.png', flush=True)
    url = os.path.join(LiCSARweb, trackID, frameID, 'metadata', 'network.png')                             # construct the URL for the metadata file
    download_data(url, outdir / 'network.png')                                                             # 




    
#%%

def get_folders_from_LiCSAR_portal(url, double_date = True):
    """On the LiCSAR S1 hub, there are folders of links to each product (e.g. for each epoch the radar intensity, or for interferograms between dates)
    This function returns a list of all these links/dates, and is designed to be used to creat paths to then download these.  
    Inputs:
        url | string | url of webpage with links on.  
        double_date | boolean | True if date/link is in interferogram form (YYYYMMDD_YYYYMMDD), False if in date form (YYYYMMDD)
    Returns:
        imdates_all | list of strings | dates/links of interst.  
    History:
        2021_03_26 | MEG | Modified from a LiCSBAS script.  
    """
    
    import requests
    from bs4 import BeautifulSoup
    import re
    
    response = requests.get(url)                                                                                # get the webpage with all the epochs on
    response.encoding = response.apparent_encoding                                                              # avoid garble?
    html_doc = response.text                                                                                    # turns the text into a giant string.  
    soup = BeautifulSoup(html_doc, "html.parser")                                                               # break the string into something that looks more like normal html
    tags = soup.find_all(href=re.compile(r"\d{8}"))                                                             # get just the links (which are epochs) from the html
    if double_date:
        character_end = 17                                                                                      # double date is of form YYYYMMDD_YYYYMMDD to need 17 characters.  
    else:
        character_end = 8                                                                                       # single date is of form YYYYMMDD so need 8 characters
    imdates_all = [tag.get("href")[0:character_end] for tag in tags]                                            # get the epochs (as a list)
    # _imdates = np.int32(np.array(imdates_all))                                                                # convert to numpy array
    # _imdates = (_imdates[(_imdates>=startdate)*(_imdates<=enddate)]).astype('str').tolist()                   # back to list?
    return imdates_all


#%%


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
    import numpy as np
    import os
    import sys
    
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



#%%


def comp_size_time(file_remote, file_local):
    """
    Compare size and time of remote and local files.
    Returns:
        0 : Size is identical and local time is new
        1 : Size is not identical
        2 : Size is identical but remote time is new
        3 : Remote not exist
    From LiCSBAS
    """
    import requests
    import datetime as dt
    import dateutil
    import os

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

#%%

def check_exist_wrapper(args):
    """
    Returns :
        0 : Local exist, complete, and new (no need to donwload)
        1 : Local incomplete (need to re-donwload)
        2 : Local old (no need to re-donwload)
        3 : Remote not exist  (can not compare, no download)
        4 : Local not exist (need to download)
    From LiCSBAS
    """
    
    import os

    i, n_data, url_data, path_data = args
    bname_data = os.path.basename(path_data)
    

    if os.path.exists(path_data):
        rc = comp_size_time(url_data, path_data)
        if rc == 1:
            print("Size of {} is not identical.".format(bname_data), flush=True)
        elif rc == 2:
            print("Newer {} available.".format(bname_data), flush=True)
        return rc
    else:
        return 4
    
