from tifffile import imread
import numpy as np
import scipy
import glob

import numpy as np
import javabridge
import bioformats
from bioformats import *
import warnings
import copy
import pickle
from tifffile import TiffWriter
import matplotlib.pyplot as plt
import matplotlib

def import_tif (file_head):
    tifs = glob.glob(str(file_head+"*.TIF"))
    img=[]
    for tif in tifs:
        img.append(imread(tif))
    return (np.array(img))



def start_bioformats():
    """Start a JVM containing Bioformat. Packed by Cell-Profiler"""
    javabridge.start_vm(class_path=bioformats.JARS,run_headless=True)
    print("Bioformats JVM started!")

def stop_bioformats():
    """Stop the JVM created by start_bioformats"""
    javabridge.kill_vm()
    print("Bioformats JVM ended!")

class MicMetadata:
    """A class for processing and storing OME-XML microscopic images meta data."""
    def _meta_init(self):
        key_values=["name","file_path","size_x","size_y","size_z","size_c","x_resolution",
        "y_resolution","pixel_size","z_resolution","pixel_type"]
        return dict.fromkeys(key_values)
    def __init__(self,image_path=None):
        self._metaData=self._meta_init()
        if image_path is not None:
            self.importMeta(image_path)

        bioformats.clear_image_reader_cache()
        javabridge._javabridge.reap()
    def importMeta(self, image_path):
        if javabridge.get_env() is None:
            start_bioformats()
        self.xml=omexml.OMEXML(get_omexml_metadata(path=image_path))
        self._metaData.update({"xml_original":self.xml.image()})
        self._metaData.update({"name":self.xml.image().Name})
        self._metaData.update({"file_path":image_path})
        self._metaData.update({"size_x":self.xml.image().Pixels.SizeX})
        self._metaData.update({"size_y":self.xml.image().Pixels.SizeY})
        self._metaData.update({"size_z":self.xml.image().Pixels.SizeZ})
        self._metaData.update({"size_c":self.xml.image().Pixels.SizeC})
        self._metaData.update({"x_resolution":self.xml.image().Pixels.PhysicalSizeX})
        self._metaData.update({"y_resolution":self.xml.image().Pixels.PhysicalSizeY})
        self._metaData.update({"z_resolution":self.xml.image().Pixels.PhysicalSizeZ})
        self._metaData.update({"DimensionOrder":self.xml.image().Pixels.DimensionOrder})

        if self.xml.image().Pixels.PhysicalSizeX!= self.xml.image().Pixels.PhysicalSizeY:
            warnings.warn("X and Y resolutions are not the same", UserWarning)
        else:
            self._metaData.update({"pixel_size":self.xml.image().Pixels.PhysicalSizeX})
        self._metaData.update({"pixel_type":self.xml.image().Pixels.PixelType})
        # print (self._metaData)
    def meta(self,key=None):
        if self._metaData is None:
            print("Use importMeta method to import metadata first")
        elif key is None:
            return (self._metaData)
        elif key not in self._metaData:
            print("provided key not available. Available keys are:")
            print(self._metaData.keys())
        else:
            return(self._metaData[key])



class MicImage(MicMetadata):
    """A class for storing Microscopic images and their metadata."""
    def __init__(self,image_path=None):
        super().__init__(image_path)
        if image_path is not None:
            self.pixels=np.zeros((self._metaData["size_z"],self._metaData["size_y"],self._metaData["size_x"],self._metaData["size_c"]),
                dtype=self.xml.image().Pixels.PixelType)
            self._metaData["DimensionOrder"] = "ZYXC"
            self._metaData["xml_original"].Pixels.DimensionOrder = "ZYXC"
            with ImageReader(image_path) as rdr:
                for channel in range(self._metaData["size_c"]):
                    for z_index in range(self._metaData["size_z"]):
                        self.pixels[z_index,:,:,channel]= rdr.read(c=channel,z=z_index,rescale=False)
                rdr.close()
            bioformats.clear_image_reader_cache()
            javabridge._javabridge.reap()
            self.sumprj=(np.sum(self.pixels,axis=0))
            self.maxprj=(np.amax(self.pixels,axis=0))
            self.meanprj=(np.mean(self.pixels,axis=0))
    def prj(self,method):
        """
        Computes and stores a projection of image along z axis.

        Parameters
        ----------
        self : MicImage class
        method : str
            the projection method to be used. "max", "sum", "mean".

        Returns
        ----------
        Store the projected image as maxprj, sumprj, or meanprj.
        """
        valid_methods={"max","sum","mean"}
        if method not in valid_methods:
            warnings.warn("Projection method is not valid, pick from following valid methods: {valid_methods}", UserWarning)
        elif method=="max":
            self.maxprj=(np.amax(self.pixels,axis=0))
        elif method=="sum":
            self.sumprj=(np.sum(self.pixels,axis=0))
        elif method=="mean":
            self.meanprj=(np.mean(self.pixels,axis=0))
    def getPRJ(self,method):
        valid_methods={"max","sum"}
        if method not in valid_methods:
            warnings.warn("Projection method is not valid, pick from following valid methods: {valid_methods}", UserWarning)
        elif method=="max":
            try: #check if maxprj exist
                return self.maxprj
            except:
                self.prj("max")
                return self.maxprj
        elif method=="sum":
            try: #check if sumprj exist
                return self.sumprj
            except: 
                self.prj("sum")
                return self.sumprj
        elif method=="sum":
            self.sumprj=(np.sum(self.pixels,axis=0))
    def view_max(self,fig_h=3,ch_names=None,title=None):
        """
        Plots the max projection of all channels along z-axis.
        This is a quick way of checking what image we are dealing with.

        Parameters
        ----------
        fig_h : int, the height of returned image
        ch_names: list, list of characters to be used naming the images.
        Defaults to "channel1",...
        title: str, the title of the image to be printed on top of the image.

        Returns
        ----------
        One row max-projected images in which each column contain the 
        image of one channel.
        """
        fig,ax = plt.subplots(1,self._metaData["size_c"],figsize=(fig_h*self._metaData["size_c"],fig_h))
        for i in range(self._metaData["size_c"]):
            ax[i].imshow(self.maxprj[...,i],norm=matplotlib.colors.Normalize())
            if ch_names is None:
                ax[i].set_title("channel"+str(i+1))
            if ch_names:
                ax[i].set_title(ch_names[i])
        if title:
            fig.suptitle(title, fontsize=16)
        # return(fig)
    def crop(self,channel,center_coord,crop_size,z_coord=None,z_size=1):
        """
        Crops a specific channel of an image with updated MetaData.
        Parameters
        ----------
        channel : int, channel to be cropped.
        center_coord: list of integers, the crop center along x-y axes.
        crop_size: int, the size of the square to be used for cropping.
        z_coord: int, the crop center along the z axis. Defaults to None,
        which crops only in x-y axes.
        z_size: The number slices around the z_coord to be included in the crop image.
        Deafults to 1, which ony returns one single slice.

        Returns
        ----------
        A new MicImage object cropped.
        """        
        x1=center_coord[0]-int(crop_size/2)
        x2=x1+crop_size
        y1=center_coord[1]-int(crop_size/2)
        y2=y1+crop_size
        img_crop=MicImage()
        img_crop._metaData={**self._metaData}
        img_crop.xml=self.xml


        if z_coord is not None and z_size>1:
            z1=z_coord-int(z_size/2)
            if z1<0:
                z1=0
            z2=z1+z_size
        if (z_coord is not None and z_size==1):
            z1=z_coord
            z2=z1+1
        if z_coord is None:
            z1=0
            z2=-1

        img_crop.pixels= self.pixels[z1:z2,x1:x2,y1:y2,channel]
        
        if img_crop.pixels.shape[0]==1:
            img_crop.pixels=np.squeeze(img_crop.pixels)
            img_crop.sumprj=np.squeeze(img_crop.pixels)
            img_crop.maxprj=np.squeeze(img_crop.pixels)
        else:
            img_crop.prj("max")
            img_crop.prj("sum")
        img_crop._metaData.update({"size_x": crop_size})
        img_crop._metaData.update({"size_x": crop_size})

        return img_crop
    def slice(self,x=None,y=None,z=None,c=None):
        """
        Creates a slice of the MicImage an updates MetaData.
        Parameters
        ----------
        x: slice class, marking the selection along x-axis.
        Defaults to slice(0,self._metaData["size_x"]) which 
        returns all the values along x-axis.
        y: slice class, marking the selection along y-axis. 
        Defaults to slice(0,self._metaData["size_y"]) which 
        returns all the values along y-axis.
        z: slice class, marking the selection along z-axis. 
        Defaults to slice(0,self._metaData["size_z"]) which 
        returns all the values along z-axis.
        c: slice class, marking the selection along c-axis. 
        Defaults to slice(0,self._metaData["size_c"]) which 
        returns all the channels.

        Returns
        ----------
        A sliced MicImage class with updated MetaData.
        """

        if x is None:
            x = slice(0,self._metaData["size_x"])
        if y is None:
            y = slice(0,self._metaData["size_y"])
        if z is None:
            z = slice(0,self._metaData["size_z"])
        if c is None:
            c = slice(0,self._metaData["size_c"])

        img_crop=MicImage()
        img_crop._metaData={**self._metaData}
        img_crop.xml=self.xml


        img_crop.pixels= copy.deepcopy(self.pixels[z,x,y,c])
        
        if img_crop.pixels.shape[0]==1:
            img_crop.pixels=np.squeeze(img_crop.pixels)
            img_crop.sumprj=np.squeeze(img_crop.pixels)
            img_crop.maxprj=np.squeeze(img_crop.pixels)
        else:
            img_crop.sumprj=(np.sum(img_crop.pixels,axis=0))
            img_crop.maxprj=(np.amax(img_crop.pixels,axis=0))
            img_crop.meanprj=(np.mean(img_crop.pixels,axis=0))

        img_crop._metaData.update({"size_z": img_crop.pixels.shape[0]})
        img_crop._metaData.update({"size_x": img_crop.pixels.shape[1]})
        img_crop._metaData.update({"size_y": img_crop.pixels.shape[2]})
        img_crop._metaData.update({"size_c": img_crop.pixels.shape[3]})

        return img_crop

    def save_ome(self,file_path):
        with TiffWriter(file_path) as tif:
            tif.save(np.moveaxis(self.pixels,[3],[1]),
                    photometric='minisblack',
                    metadata={'axes': 'ZCYX'})# need to handle meta data batter.

    def pickle(self,file_path):
        with open(file_path, "wb") as f:
            pickle.dump(self, f)


def tiff_imp(file_path,zxy_ch):
    '''
    Import a ome.tif file and returns it as an numpy ndarray.
    
    It move the first axis to the end to match skimage dimension recommendations.
    
    Prameters:
    file_path: Location of the ome.tif file.
    zxy_ch: NEED explanation

    Returns:
    ndarray: Pixel values of ome.tif file.

    '''

    from skimage import io
    import numpy as np
    img = io.imread(file_path)
    #In my images the channel is the first dimension
    #skikit recommends to have the following order for fast computation (Z,X,Y,ch)
    return np.moveaxis(img,zxy_ch,[0,1,2,3])

def dv_import(file_path,ij):
    '''
    Imports a DeltaVision file using FIJI. This requires an imagej Java VM to be open.

    Prameters:
    file_path: location of .dv file
    ij: imageJ instance

    Returns:
    Image as an numpy array. The axes order is z,x,y,ch.
    '''

    import numpy as np

    img=ij.io().open(file_path)
    as_np=ij.py.from_java(img)
    np_reorder=np.moveaxis(as_np,0,-1)

    return np_reorder

def zstack_import(file_path,ij,file_type):
    '''
    Imports a DeltaVision file using FIJI. This requires an imagej Java VM to be open.

    Prameters:
    file_path: location of .dv file
    ij: imageJ instance

    Returns:
    Image as an numpy array. The axes order is z,x,y,ch.
    '''

    import numpy as np

    img=ij.io().open(file_path)
    # img=ij.openImage(str(file_path))
    as_np=ij.py.from_java(img)
    if file_type=="dv":
        np_reorder=np.moveaxis(as_np,0,-1)
    elif file_type=="flex":
        np_reorder=np.moveaxis(as_np,1,-1)
    # img.close()
    return np_reorder


