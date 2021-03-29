# LiCSAR Web Tools

A small collection of Python3 scripts to download data from the COMET-LiCS Sentinel-1 InSAR portal.  The licsar_web_tools_example.py shows two examples, but usage can be summaried below:

<pre><code>
frameID = '014A_04939_131313'                                       # As defined on the LiCSAR portal
date_start = 20190901                                               # YYYYMMDD
date_end   = 20200101                                               # YYYYMMDD
download_metadata = True                                            # the network of baselines, DEM etc.  
between_epoch_files = ['geo.unw.png', 'geo.cc.png']                 # Possible files: 'geo.cc.png', 'geo.cc.tif', 'geo.diff.png', 'geo.diff_pha.tif', 'geo_diff_unfiltered.png', 'geo_diff_unfiltered_pha.tif', 'geo.unw.png', 'geo.unw.tif'
epoch_files = ['geo.mli.png','ztd.jpg']                             # Possible files: 'geo.mli.png', 'geo.mli.tif', 'sltd.geo.tif', 'ztd.geo.tif', 'ztd.jpg'
n_para = 4                                                          # Parallelisation.  The number of cores is a good starting point.    


download_LiCSAR_portal_data(frameID, date_start, date_end, download_metadata, epoch_files, between_epoch_files, n_para)

</code></pre>


This has only been tested on Linux (Ubuntu 18.04), so please get in touch if you encounter any bugs.  Several of the functions and ideas used in this software were  were taken from LiCSBAS: <https://github.com/yumorishita/LiCSBAS>
