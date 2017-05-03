from PyQt4.QtCore import *
from qgis.core import *
import gdal
import numpy as np

#import doViewshed
# this is circular import !! cannot do import dem_chunk

#modes to combine outputs TO IMPLEMENT !!!!!!!!!!!!!!!!!!!!!!!!
SINGLE = 0
ADD = 1
MIN = 2
MAX = 3
# ------------------------------------------

'''
This class handles input and output of raster data.
'''

class Raster:

    
    def __init__(self, raster, crs=None):
	
		
        gdal_raster=gdal.Open(raster)
        
        
        # This is problematic : the idea is that the thing might work even with 
        # rasters not loaded in QGIS. But for those already loaded, it will override gdal crs. 
        # format has to be same (text string) ??
       
        self.crs = crs if crs else gdal_raster.GetProjection()
        
        self.rst = gdal_raster #for speed, keep open raster ?
                        
        self.mode= 0 #buffer mode is none
        # size is y first, like numpy
        self.size = (gdal_raster.RasterYSize, gdal_raster.RasterXSize)
           
    
        #adfGeoTransform[0] /* top left x */
        #adfGeoTransform[1] /* w-e pixel resolution */
        #adfGeoTransform[2] /* rotation, 0 if image is "north up" */
        #adfGeoTransform[3] /* top left y */
        #adfGeoTransform[4] /* rotation, 0 if image is "north up" */
        #adfGeoTransform[5] /* n-s pixel resolution */

        gt=gdal_raster.GetGeoTransform()
        self.pix = gt[1] 
                
        raster_x_min = gt[0]
        raster_y_max = gt[3] # it's top left y, so maximum!

        raster_y_min = raster_y_max - self.size[0] * self.pix
        raster_x_max = raster_x_min + self.size[1] * self.pix

        
        self.extent = [raster_x_min, raster_y_min, 
                       raster_x_max, raster_y_max]


    """
    Determine the method for combining results (summing up or otherwise).
    By default, the mode is 0: no summing up
    Buffer size (self.result) is the same size as the entire array
    [ should be done in chunks for very large arrays - to implement ...]
    """

    def set_buffer_mode(self, mode):

        self.mode = mode
        self.result = np.zeros(self.size)

        
        

    """
    This is the largest window, used for all analyses.
    Smaller windows are slices of this one.

    [ theoretically,window should be a subclass,
    but we can have only one window at a time,
    and it should be as fast as possible] 
    """
    def set_master_window (self, radius_pix):           

  
        full_size = radius_pix *2 +1
        self.window = np.zeros((full_size, full_size))

      
    def get_diameter_earth (self):

        crs= self.crs		
    
        start = crs.find("SPHEROID") + len("SPHEROID[")
        end = crs.find("]],", start) + 1
        tmp = crs[start:end].split(",")

        try:
                semiMajor = float(tmp[1])
                if not 6000000 < semiMajor < 7000000:
                        semiMajor = 6378137
        except:
                semiMajor = 6378137

        try:
                flattening = float(tmp[2])
                if not 296 < flattening < 301:
                        flattening = 298.257223563
        except:
                flattening = 298.257223563

        semiMinor = semiMajor - semiMajor / flattening
        
        return semiMajor + semiMinor

 
    """
    TODO : emit a warning ?!
    """
    def set_mask (self, mask):

        if mask.shape == self.window.shape: self.mask = mask
        else : pass


    """
    Extract a quadrangular window from the raster file.
    Observer point (x,y) is always in the centre.

    Upon opening a window, all parameters regarding its size and position are
    registered in the Raster class instance - and reused for writing results
    """
    def open_window (self, x, y, radius_pix):

        rx = radius_pix
        #to place smaller windows inside the master window
        diff_x = self.window.shape[1] - (rx *2 +1)
        diff_y =  self.window.shape[0] - (rx *2 +1)

        

        if x <= rx:  #cropping from the front
            x_offset =0
            x_offset_dist_mx = (rx - x) + diff_x
        else:               #cropping from the back
            x_offset = x - rx
            x_offset_dist_mx= 0

                           
        x_offset2 = min(x + rx +1, self.size[1]) #could be enormus radius, so check both ends always
        
        if y <= rx:
            y_offset =0
            y_offset_dist_mx = (rx - y) + diff_y
        else:
            y_offset = y - rx
            y_offset_dist_mx= 0

        y_offset2 = min(y + rx + 1, self.size[0] )

        window_size_y = y_offset2 - y_offset
        window_size_x = x_offset2 - x_offset

        self.window_slice = np.s_[y_offset : y_offset + window_size_y,
                                  x_offset : x_offset + window_size_x ]


        self.inside_window_slice = np.s_[ y_offset_dist_mx : y_offset_dist_mx +  window_size_y,
                                          x_offset_dist_mx : x_offset_dist_mx + window_size_x]

        self.gdal_slice = [x_offset, y_offset, window_size_x, window_size_y]

 
        
        self.window [:]=0 # as a precaution, not needed ??

        self.window[ self.inside_window_slice] = \
                         self.rst.ReadAsArray(*self.gdal_slice ).astype(float)


##        self.offset = (x_offset, y_offset)
##        self.win_offset= (x_offset_dist_mx, y_offset_dist_mx)
##        self.win_size = (window_size_x, window_size_y)




    def add_result(self, in_array):

        if self.mode == 0: self.result = in_array

        elif self.mode == 1:
            
            self.result [self.window_slice] += in_array [self.inside_window_slice]
            


        elif self.mode == 2:
            np.where(self.result [self.window_slice] < in_array [self.inside_window_slice],
                     in_array [self.inside_window_slice], self.result [self.window_slice])

        elif self.mode == 3:
            print  '-------- TODO ! ---------'
        
            

    """
    TODO : a trick to work on very large arrays [mode = cumulative_lage]
    e.g. read a window from a gdal raster, sum, write back
       
    """
    def write_result(self, filePath,
                     fill = np.nan, no_data = np.nan,
                     dataFormat = gdal.GDT_Float32):
        
            
        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create(filePath, self.size[1], self.size[0], 1, dataFormat)
        ds.SetProjection(self.crs)
        ds.SetGeoTransform(self.rst.GetGeoTransform())


        ds.GetRasterBand(1).Fill(fill)
        ds.GetRasterBand(1).SetNoDataValue(no_data)

##        #cumulative : same size as original array
##        if self.mode == 'cumulative':  ds.GetRasterBand(1).WriteArray(self.result)
##
##        else:
        
        #all modes > 0 operate on a copy of the raster
        if self.mode:
            
            ds.GetRasterBand(1).WriteArray(self.result)
        else:
            #for writing it takes only x and y offset (1st 2 values of self.gdal_slice) 
            ds.GetRasterBand(1).WriteArray(self.result[ self.inside_window_slice ],
                                           *self.gdal_slice[:2] )
            #self.offset[0], self.offset[1])
        ds = None

	  

    
        

