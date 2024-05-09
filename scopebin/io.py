from tifffile import imread
import numpy as np
import scipy
import glob

def import_tif (file_head):
    tifs = glob.glob(str(file_head+"*.TIF"))
    img=[]
    for tif in tifs:
        img.append(imread(tif))
    return (np.array(img))