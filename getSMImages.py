'''
        ***  DOWNLOAD SOLAR IMAGES FROM SOLARMONITOR ***
'''

__author__ = " Gabriel García García "

__email__  = " gabgarar@gmail.com "


from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange
import urllib.request
import os
import gzip
import shutil

def getURLFromSM(inst,wv,soon,late):
     """
        @ params
            * inst : instrument 
            * wv : length wave 
            * soon  : date more recent
            * late  : date more late
        
        @ return
            * ret : list with urls from SM in that range of time
        
    """
    solmon_pattern = (
                             'http://solarmonitor.org/data/'
                             '%Y/%m/%d/fits/{instrument}/'
                             '{instrument}_{wave:05d}_fd_%Y%m%d_%H%M%S.fts.gz'
                       );
    solmon = Scraper(solmon_pattern, instrument = inst, wave = wv);
    
    timerange = TimeRange(late,soon); 
    list_urls = solmon.filelist(timerange); 
    return list_urls; 


def downloadFITSfrom(urls,dir_cont):

    """
        @ params
            * urls : list with urls from SM
            * dir_cont : name of the dir where we put the images
        
        @ return
            * ret : all images decompressed in dir_cont
        
    """
    print(" The container folder is : " + dir_cont + "\n") ;
    try:
      os.stat(dir_cont) ; 
    except:
      os.mkdir(dir_cont) ; 
    
    # The first one is download the archives into dir_cont 
    
    print("Downloading ..." + "\n");
    
    for elem in urls:
        dir_out = elem[elem.index(dir_cont):] ; 
        
        print("\t" + "Path : " + dir_out + "\n"); 
            
        archivo_tmp, header = urllib.request.urlretrieve(elem)
        with open(dir_out, 'wb') as archivo:
           with open(archivo_tmp, 'rb') as tmp:
                archivo.write(tmp.read())
      
    # Second the archives are decompressed
    
    print("Decompressing ..." + "\n")
    for path_fits in os.listdir(dir_cont):
            if ( path_fits[-3:] == '.gz'):
                path     = dir_cont + "/" + path_fits ; 
                print("\t" + "The archieve will be descompressed as :" + path[:-3] + "\n"); 
                with gzip.open(path, 'rb') as f_in:
                    with open(path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out);
                        f_in.close(); 
                        f_out.close();
                        os.remove(path); 
                        
                        
        