'''
        ***  GENERATE SYNOPTIC SOLAR MAP ***
'''

__author__ = " Gabriel García García "

__email__  = " gabgarar@gmail.com "


# Importing libraries

import os
import numpy as np
import matplotlib.pyplot as plt
import cv2

import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

import sunpy.physics
import sunpy.cm
from sunpy.coordinates import frames
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate,differential_rotate

import squircle
import sys
import argparse




def getPointToPixelWorld(alf, bet, aia_map):
    
    """
        @ params
            * alf : in degrees
            * bet : in degrees
            * aia_map : structure individual map
        
        @ return
            * ret : Convert a point in HelioGraphicsStonyHurst to pixels
        
    """
    point = SkyCoord(alf * u.deg,bet * u.deg,frame=frames.HeliographicStonyhurst);
    point = point.transform_to(aia_map.coordinate_frame)
    return aia_map.world_to_pixel(point);

  
def cropImageEdge(aia_map):
    
    """
        @ params
            * aia_map : structure individual map
        
        @ return
            * ret : Returns an image cropped by the edge of the sun
        
    """
    aba     = getPointToPixelWorld(-90,0, aia_map)
    arr     = getPointToPixelWorld(90,0, aia_map)
    der     = getPointToPixelWorld(0,90, aia_map)
    izq     = getPointToPixelWorld(0,-90, aia_map)
    point   = izq.y, arr.x;
    limites = point[0].astype(int)/u.pix,2*point[0].astype(int)/u.pix,point[1].astype(int)/u.pix,2*point[1].astype(int)/u.pix

    return aia_map.data[limites[0]:limites[2],limites[0]:limites[2]];



def cropImageNmiddle(img_exp,aux,N):

    """
        @ params
            * img_exp : Expanded array image
            * aux : How many pixels you want to fill between two images to reduce the posible error.
            * N : 
        
        @ return
            * ret : Returns an expanded and cropped image in the middle
        
    """            
    pto = N * ((img_exp.shape[0] / 180) * 13.19)+ aux;
    pto = int(pto);

    ini = ( img_exp.shape[1] - pto ) / 2; 
    
    return img_exp[0:img_exp.shape[0],int(ini - 1): int((ini - 1) + pto)];



def mappingTheImage(path_dir,path_fits):
    
    """
        @ params
            * path_dir : Directory where the images are located
            * path_fits : Actual image
        
        @ return
            * ret : map from actual image
        
    """ 
    hdu_list_princ   = fits.open(path_dir + '/' + path_fits);
    header_princ     = hdu_list_princ[0].header;
    img_data_princ   = hdu_list_princ[0].data;

    return sunpy.map.Map((img_data_princ,header_princ)) ;## Tener cuidado ya que no se si usar img_expand o img_data.


def mergeSolarImagesFromDir(path_dir,numDays,path_out,aux,obser):
    
    '''
        @ params
            * path_dir : Directory where the images are located
            * numDays : How many days you want to represent into your synoptic map
            * path_out : Where you want to save the synoptic map
            * aux : How many pixels you want to fill between two images to reduce the posible error.
            * obser : What instrument we are using to get the cmap to give pseudocolor

        @ return
            * ret : Get a synoptic map from numDays days from your images.

    '''

    
    cont = 0; 
    cont_temp = 0; 

    for path_fits in os.listdir(path_dir):

        if (cont_temp > numDays):
           break;
        

        # map the image
        aia_map    =  mappingTheImage(path_dir,path_fits)

    
        
        if (os.listdir(path_dir)[0] == path_fits): 
            N = 1;
            date = aia_map.date; 
            print("\t\t[",cont,", 0 ] Base image loaded.");
    
        else:
            N = (aia_map.date - date).value 
            cont_temp += N; 

        if ( N >= 0):
            
            # We cut by the edge of the sun
            crp = cropImageEdge(aia_map); 
            
            # Expand the image
            img_exp = squircle.to_square(crp, method="fgs");
            
            # Crop the image through the center
            img_pq  = cropImageNmiddle(img_exp,aux,N); 
            

            # Put he first image as base image
            if (os.listdir(path_dir)[0] == path_fits):
                final       = img_pq.copy() ; 
                tamX, tamY  = img_pq.shape;

            else :
                img_pq  = img_pq[0:tamX,0:tamY];
                print("\t\t[",cont,", 1 ] Concatenating images.");
                final    = np.concatenate((img_pq,final),axis=1);
                print("\t\t[",cont,", 2 ] Composite image with success.");

            date = aia_map.date ;
            cont += 1
            
    print("::::::::::::::::::::::::::::::::::::::::::::");
    print("Image will be saved in : ", path_out);
    print("Saving...");
    plt.imsave(path_out,final,cmap=plt.get_cmap(obser));  
    
    