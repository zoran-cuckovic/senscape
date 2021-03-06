# -*- coding: utf-8 -*-

"""
/***************************************************************************
ViewshedAnalysis
A QGIS plugin
begin : 2013-05-22
copyright : (C) 2013 by Zoran Čučković
email : /
***************************************************************************/

/***************************************************************************
* *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation version 2 of the License, or *
* any later version. *
* *
***************************************************************************/
"""

from __future__ import division

"""
BUGS
- progress bar not working: nicer [NOT VITAL]
- finding highest position : should be circular [NOT VITAL]
- Limitations in cell size eg. : SAGA assumes that raster layers have the same cell size in the X and Y axis. If you are working with a layer with different values for horizontal and vertical cell size, you might get unexpected results. In this case, a warning will be added to the processing log,
  indicating that an input layer might not be suitable to be processed by SAGA
"""

from PyQt4.QtCore import *
#from PyQt4.QtGui import *
from qgis.core import *
#from qgis.utils import * #progress bar
import os
#from osgeo import osr, gdal, ogr

import time
#_____Testing_____________
#from cProfile import Profile

import numpy as np
from math import sqrt

from Points import Points
from Raster import Raster

#from processing.gui.SilentProgress import SilentProgress NOT WORKING

BINARY_VIEWSHED = 0
INVISIBILITY_DEPTH = 1
HORIZON = 2
HORIZON_FULL = 3
ANGULAR_SIZE = 4
# not used: separate function!
INTERVISIBILITY = 5



def dist(x1,y1,x2,y2, estimation=False):
    if not estimation: r=sqrt(pow((x1-x2),2) + pow((y1-y2),2))
    else: # error = cca 1% - NOT USED!
        rt= 1.4142135623730950488016887242097
        r = (abs(x1-x2) * rt + abs(y1-y2) * rt) / 2

    return r

"""
Return a map of distances from the central pixel.
Attention: these are pixel distances, not geographical !
(to convert to geographical distances: multiply with pixel size)
"""
def distance_matrix ( radius_pix, squared=False):
          
    full_window_size = radius_pix *2 + 1
    
    temp_x= ((np.arange(full_window_size) - radius_pix) ) **2
    temp_y= ((np.arange(full_window_size) - radius_pix) ) **2

    if not squared:
        return np.sqrt(temp_x[:,None] + temp_y[None,:])
    # squared values
    else: return temp_x[:,None] + temp_y[None,:]

"""
Model vertical drop from a plane to spherical surface (of the Earth)
Note that it has to be multiplied with pixel size to get usable values

"""
def curvature_matrix(radius_pix, diameter_earth, pixel_size, refraction):

    dist_squared = distance_matrix(radius_pix, squared=True)
    # all distances are in pixels in doViewshed module !!
    # formula is  squared distance / diam_earth 
    # need to divide all with pixel size (squared !!)
    D = diameter_earth / (pixel_size **2)
        
    return (dist_squared / D) * (1 - refraction) 

def error_matrix(radius, size_factor=1):

    """
    Create a set of lines of sight which can be reused for all calculations. 
    Each line (departing from the observer point) has its target and a set of pixels it passes through.
    Only 1/8th of full radius is enough : the rest can be copied/mirrored. 
    """

    if size_factor == 0: size_factor = 1 #0 is for non-interpolated algo...
    radius_large = radius  * size_factor  
                                                
    mx_index= np.zeros((radius_large +1 , radius, 2)).astype(int)
    mx_err = np.zeros((radius_large +1 , radius))
    mx_mask = np.zeros(mx_err.shape).astype(bool)

    min_err = {}

    j=0 #keep 0 line empty

    for m in range (0, radius_large+1 ): # 45 deg line is added (+1) 

        x_f, y_f = radius, radius #x0,y0

        #dy = x; dx = y 
        dy,dx= m, radius_large #SWAPPED x and y! MESSY

    
        #x and y = delta x and y but y is steep!
        #fist line is min y then it ascends till 45°

        D=0
        for i in xrange (0, radius ):   #restrict iteration to actual radius!     
            x_f += 1
            if 2* (D + dy) < dx:
                D += dy # y_f remains
            else:
                y_f += 1
                D += dy - dx
                           
            #reverse x,y for data array!
            yx= (y_f,x_f)
            mx_index[j,i,0:2]=yx
                           
            if D: e=D/dx; err=abs(e)
            else: e, err = 0,0

            mx_err[j,i]=e
          # keep pixel dictionary to sort out best pixels
            try:
                err_old = min_err[yx][0] 
                if err < err_old: min_err[yx]=[err,j,i]
            except:
                min_err[yx]=[err,j,i]
   
        j+=1
    

    #check-out minimum errors
    # numpy style would be np.argmin.at!
    for key in min_err:
        ix=min_err[key][1:3]
        er = min_err[key][0]
        mx_mask[ix[0], ix[1]]= 1

    mx_err_dir = np.where(mx_err > 0, 1, -1)
    mx_err_dir[mx_err == 0]=0 #should use some multiple criteria in where... 

    #take the best pixels  
    #cannot simply use indices as pairs [[x,y], [...]]- np thing...
    #cannot use mx : has a lot of duplicate indices

   
    mx_err_index = mx_index [:,:, 0] + mx_err_dir
                                # we do not need negative errors any more
    return mx_index, mx_err_index, np.absolute(mx_err), mx_mask


"""
    Single point viewshed function: works on a number of precalculated matrices: for speed.
    Takes prepared errors and indices of best pixels (with least offset from line of sight)
    Cannot be much simplified (?) - without loosing speed...
    Note that only 1/8 of the entire analysed zone is enough for interpolation etc,
    the rest is filled by flipping arrays.
    
"""
    
def viewshed_raster (option, data,                  
              error_matrix, error_mask, indices, indices_interpolation,                  
              target_matrix=None, distance_matrix=None, interpolate = True):


    #np.take is  much more efficient than using "fancy" indexing (stack overflow says ...)
    #... but it flattens arrays ! TODO...
    
    center = int(data.shape[0] / 2) #ugly, but passing is as argument is even uglier...
    # could be ... data.window_center, but data should be Raster object ...
    
    mx_vis = np.zeros(data.shape)#if np.ones then the centre gets True val (for Binary)
   
    mx = indices[: ,:, 1]; my = indices[: ,:, 0]

    me_x, me_y = mx, indices_interpolation #me_x = mx !

    # flipping ang flopping data so that it fits indices matrix
    # (rather than producing multiple indices matrices :) )
    views = [ np.s_[:], np.s_[ :, ::-1],
              np.s_[ ::-1, :], np.s_[ ::-1, ::-1] ]
       

    for steep in [False, True]: #- initially it's steep 0, 0

        if steep: #swap x and y
            
            me_x, me_y = me_y , me_x 
            mx, my = my,mx 

   
        
        for view in views:                

            view_d = data[view]
            view_o = mx_vis[view]              
           
            interp = view_d[mx,my]
           

            if interpolate: 
                interp += (view_d[me_x, me_y] - interp) * error_matrix
           
                           
            # do it here so we can subsitute target below!
            test_val = np.maximum.accumulate(interp, axis=1)
            
            
            if target_matrix is not None:

                view_tg = target_matrix[view]
                # substitute target matrix, but only after test matrix is calculated!
                interp = view_tg[mx,my] 

                if interpolate: 
                
                # could be done only on "good pixels", those with minimal error !!
                # use mask on mx, my etc  - test for speed ...
                    interp += (view_tg[me_x, me_y] - interp) * error_matrix
     
            # non-interpolated, normal                  
           # v = data[mx,my] == np.maximum.accumulate(data[mx,my], axis=1)

            if option == BINARY_VIEWSHED :  v = interp >= test_val
           
            elif option == INVISIBILITY_DEPTH : v = interp - test_val

            elif option == ANGULAR_SIZE:
                v = np.zeros(interp.shape)#This is INEFFICIENT           
                v[:, 1:] = np.diff(test_val)

                
            elif option == HORIZON: # = true or last horizon

                v = interp >= test_val               
                                
                                           
                #select last pixel (find max value in a reversed array (last axis!)
                #argmax stops at first occurence
                #indices have to be re-reversed :)
                #gives a flat array of 1 index for each row (axis=1)
                rev_max = center - np.argmax(v[:, ::-1], axis=1) -1
                               
                v[:] = False

                #radius = row n° for fancy index (should be some nicer way...               

                v[ np.arange(center +1), rev_max.flat ] = True

 

            elif option == HORIZON_FULL:
                
                v = interp >= test_val     
                #diff always returns one place less than full array! 
                v[:, :-1] = np.diff((interp >= test_val).astype(int)) *-1
                           
                # 1: make = True/False array;
                # 2: turn to integers because boolean operations give True False only
                # 3 diff = x(i)- x(i-1) = -1 or 1 for breaks
                # * -1 so that good breaks become +1, instead of -1
                                
                v[v == -1]=0 #delete the beginning of  visible areas
                            
            #mx_vis [mx[mask], my[mask]]= v[mask] #np.absolute(mx_err[mask]) for errors

            #view of mx_vis, NOT A COPY!
            view_o [mx[error_mask], my[error_mask]]= v[error_mask]
            

    # handling some details ...
    if option== BINARY_VIEWSHED: mx_vis[center,center]=True

    elif option== INVISIBILITY_DEPTH:
            # = OPTION ANGLE ZA ALGO : angle diff , angle incidence ???
                
        mx_vis *= distance_matrix 
        # assign target height to the centre (not observer height !)
        # = first neigbour that is always visible :)

        mx_vis[center,center]=mx_vis[center,center+1]

    
# ---------------------------------------        
    return mx_vis


def rasterised_line (x,y, x2, y2, interpolation = True):

    """
    Calculate intervisibilty lines from the observer point (always in the centre of the matrix)
    to target point (x2, y2).
    Has it's proper algorithm in order to avoid inaccuracies of the usual viewshed approach.
    
    """

  
    
    dx = abs(x2 - x); dy = abs(y2 - y)
    steep = (dy > dx)
    #direction of movement : plus or minus in the coord. system
    sx = 1 if (x2 - x) > 0 else -1
    sy = 1 if (y2 - y) > 0 else -1

    if steep: # if the line is steep: change x and y
        #x,y = y,x they are the same !!
        
        dx,dy = dy,dx
        sx,sy = sy,sx

    D = 0
  
    #for interpolation
   # slope = dy / dx *sx *sy #!!
   #the operators for negative quadrants (we do need dx, dy as absolute to verify steepness, to search distances...)

    dx_short = dx-1 # to leave out the last pixel (target)
    
    #store indices 1) los, 2) neighbours, 3) error
    mx_line = np.zeros((dx_short, 2), dtype=int)
    
    if interpolation:
        mx_neighbours = np.zeros((dx_short,2), dtype=int)
        mx_err = np.zeros((dx_short))
        msk = np.ones((dx_short),dtype=bool)


    for i in xrange (0, dx_short): 

    # ---- Bresenham's algorithm (understandable variant)
    # http://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html       
        x += sx
        if 2* (D + dy) < dx:
            D += dy # y remains
        else:
            y += sy
            D += dy - dx
                       
        #unfortunately, np allows for negative values...
        # when x or y are too large, the break is on try-except below
        
  
        # --------- coordinates not unpacked ! ------------
        
        mx_line[i, :] = [y, x] if not steep else [x, y] 
      
        if interpolation:
 
            if D:
                sD = -1 if D < 0 else 1
                interp = y + sy *sD

                if steep:
                    mx_neighbours[i, :] = x, interp
                else:
                    mx_neighbours[i, :] = interp, x

                mx_err [i]=D 
  
            else:   msk[i]=False

    if interpolation:
        #give zero-error points themselves as neighbours
        # NB. this is not needed because error is 0; the result *= 0,
        # but it's more clear this way, and will eliminate the possibility of stepping out of matrix
        mx_neighbours[~msk, :]= mx_line[~msk, :]
        
        #error is actually D / dx !
        mx_err[msk]/= dx # zero values will give nans on division!
        return mx_line, mx_neighbours, abs(mx_err) 
    else: return mx_line
            

"""
Calculates intervisibility lines...

"""

def intervisibility (points_class, targets_class, raster_class,
                      curvature=0, refraction=0, interpolate = True):


    ###########################"
    edges = {}

    # RasterPath= str(QgsMapLayerRegistry.instance().mapLayer(Raster_layer).dataProvider().dataSourceUri())
    points, targets = points_class.pt, targets_class.pt
    
    radius_pix = int(points_class.max_radius)
    
    raster_class.set_master_window(radius_pix)
    
    mx_dist = distance_matrix(radius_pix)
         
    if curvature:
        #everything is in pixels here, so the earth diam
        #has to accord!
        mx_curv = curvature_matrix(radius_pix,
                    raster_class.get_diameter_earth(),
                                   raster_class.pix,
                                    refraction )
        
    else: mx_curv = 0
       
    for id1 in points :

        try: tg = points[id1]["targets"]
        except: continue
        
        x,y= points[id1]["x_pix"],  points[id1]["y_pix"]
        z= points[id1]["z"]
        
        # radius is in pixels !
        r=  points[id1]["radius"]
        #r_pix= int (r)

     
        raster_class.open_window (x, y, radius_pix)

        data=raster_class.window
        
        z_abs = z + data [radius_pix,radius_pix]
        
        # level all according to observer
        data -=  z_abs + mx_curv 
    
        data /= mx_dist #all one line = (data -z - mxcurv) /mx_dist

        for id2 in tg:
            #adjust for local raster (diff x)
            x2 = radius_pix + (targets[id2]["x_pix"] - x)
            y2 = radius_pix + (targets[id2]["y_pix"] - y)

            angle_targ = data[y2,x2] 

            d = mx_dist[y2,x2]
            
            z_targ = targets[id2]["z_targ"]
            
            mx_line, mx_neighbours, mx_err = rasterised_line (
                                    radius_pix, radius_pix, x2, y2,
                                    interpolation= interpolate)

            l_x, l_y = mx_line[:,1], mx_line[:,0]

            angles = data[l_y,l_x]            
            
            if interpolate:
                n_x, n_y = mx_neighbours[:,1], mx_neighbours[:,0]
                angles +=  (data[n_y, n_x] - angles) * mx_err

            # bare terrain!
            depth = (angle_targ - np.max(angles)) * d
            
            # correct with target height only for invisible terrain,
            # otherwise it adds to relative pixel height
            edges[id1, id2] = z_targ if depth >= 0 else depth + z_targ

            print id1, id2, depth
                          
    return edges  
    
"""
Opens a DEM and produces viewsheds from a set of points. Algorithms for raster ouput are all
based on a same approach, explained at zoran-cuckovic.com/landscape-analysis/visibility/.
Intervisibility calculation, however, is on point to point basis and has its own algorithm.
"""


def Viewshed (points_class, raster_class, 
          output_options,
          Target_points=None,
          curvature=0, refraction=0, algorithm = 1,
          progress = None): 

          

    ##### TESTING #########
    start = time.clock(); start_etape=start
    test_rpt= "Start: " + str(start)


    #prof=Profile()

                   
    #prof.enable()

    out_files=[]
    rpt=[]

    # RasterPath= str(QgsMapLayerRegistry.instance().mapLayer(Raster_layer).dataProvider().dataSourceUri())
    points = points_class.pt

    # this is maximum radius in *pixel distance*
    # should be explicit ( .max_radius(pixel=True) ) ...
    radius_float = points_class.max_radius
    radius_pix = int(radius_float)

    #for speed and convenience, use maximum sized window for all analyses
    #this is not clear! should set using entire size, not radius !!
    raster_class.set_master_window(radius_pix)
    
    #pre_calculate distance matrix

    mx_dist = distance_matrix(radius_pix)
         
    if curvature:
        #everything is in pixels here, so the earth diam
        #has to be divided with pixel size!
        mx_curv = curvature_matrix(radius_pix,
                    raster_class.get_diameter_earth(),
                    raster_class.pix,
                               refraction)

       
    else: mx_curv = 0
    
    mx_indices, mx_err_indices, mx_err, mx_mask = error_matrix(radius_pix, algorithm)

    cnt = 0
    for id1 in points :

         # x,y, EMPTY_z, x_geog, y_geog = points[id1] #unpack all = RETURNS STRINGS ??
        
        x,y= points[id1]["x_pix"],  points[id1]["y_pix"]
        z= points[id1]["z"]; z_targ= points[id1]["z_targ"]
        
        # radius is in pixel units !
        r=  points[id1]["radius"]  

        # all matrices are calculated for the maximum radius !!        
        

##        if diff: np_slice = np.s_[diff : radius_pix - diff,
##                                  diff : radius_pix - diff]
##
##        else: np_slice = np.s_[:]
##        
        
        # background value for analysed matrix
        # it gives observer value, as that point is never tested

        # should not take too much time,
        # compared to (faster) data [y, x] = init_val (test?)
##        if output_options == BINARY_viewsghed: init_val = 1
##        elif output_options == INVISIBILITY_DEPTH: init_val = z_targ
##        else: init_val = 0

    
        # this routine is also asigning offsets of data window - to the Raster class!
        raster_class.open_window (x, y, radius_pix,
                                  initial_value = 0,
                                  pad = algorithm > 0)

        data=raster_class.window


        #used for calculating extracted area, for the report, and for horizon crop
        # ... actually should be a function of raster class: get_area()
        s_y, s_x = raster_class.inside_window_slice 
        
        z_abs = z + data [radius_pix,radius_pix]
        
        # level all according to observer
        data -= mx_curv + z_abs

        data /= mx_dist #all one line = (data -z - mxcurv) /mx_dist
            
        if z_targ : mx_target = data + z_targ / mx_dist 
                        # it's ( data + z_target ) / dist,
                        # but / dist is already made above
        else: mx_target=None
        
        # Horizon is the last visible zone beor the edge of the window
        # need to remove data from corners (if a circular analysis is required !)
        if output_options in [HORIZON, HORIZON_FULL]:
            data [mx_dist >= r + 2] = np.min(data)

               
        matrix_vis = viewshed_raster (output_options,  data,
                                      mx_err , mx_mask,
                                      mx_indices, mx_err_indices,
                                      target_matrix=mx_target,
                                      distance_matrix=mx_dist,
                                      interpolate = algorithm > 0)


        if output_options in [HORIZON, HORIZON_FULL]:

            #crop borders 
            #cannot solve the problem when the analysis window is larger
            #than DEM extents - the outside values are
            # forcing fake horizons on borders...

            msk = np.ones(data.shape).astype(bool)
                        
            msk[s_y[0] +1 : s_y[1] -1, s_x[0] +1 : s_x[1] -1]= False

            matrix_vis[msk]=0

        # ----------- MASKING ------------
        # TODO : complex masking ......
        mask_circ = mx_dist [:] > r


        # we need zero background to sum up,
        # but not for other options
     
        fill = 0 if raster_class.mode == 1 else np.nan
        matrix_vis[mask_circ]=fill
        
        #TODO: make some kind of general report system !
        # algo 0 is fast, so skip to save some time (??)
        if algorithm > 0 :
            
            sl = np.s_[slice(*s_y), slice(*s_x)]
            #careful with areas outside raster ! 
            view_m = matrix_vis[sl]

            if output_options == INVISIBILITY_DEPTH:
               c= np.count_nonzero(view_m >= 0) 

            else:  c = np.count_nonzero(view_m)
                 
            # Count values outside mask (mask is True on the outside!)
            crop =np.count_nonzero(mask_circ[sl])
            # Here, nans are in the outside (which are not zero!)
            if raster_class.mode != 1: c -= crop
            # this is unmasked: sunbtract masked out areas!
            rpt.append([id1, c , view_m.size-crop] )

            

            

        # do writing inside loop if required
        if raster_class.mode <= 0:
            # None will force registered file path
            if raster_class.mode == -1: f_name = None
            else: f_name = points[id1]["file"]
                        
            raster_class.write_result(file_name=f_name,
                                      in_array=matrix_vis)
            out_files.append(f_name)            
        else: raster_class.add_to_buffer (matrix_vis)



        cnt += 1
        # NOT WORKING, QGIS FREEZES
        try: progress.setPercentage( (cnt/points_class.count) *100)
        except: print "not working"
    """
    #TESTING #################
    prof.disable()
    test_rpt += "\n Total time: " + str (time.clock()- start)
    QMessageBox.information(None, "Timing report:", str(test_rpt))

    import pstats, StringIO
    s = StringIO.StringIO()
   
    ps = pstats.Stats(prof, stream=s).sort_stats('cumulative')
    ps.print_stats()
    # print s.getvalue()
    

    test_rpt += "\n Total time: " + str (time.clock()- start)
    """
    
    if raster_class.mode > 0:
        raster_class.write_result()#write using .output property
        out_files.append(raster_class.output)
    
    # TODO : error signalling (catch from raster_class...)
    return out_files, rpt



def write_intervisibility_line (file_name, data_list, coordinate_ref_system, use_pix_coords=False):

    #QMessageBox.information(None, "Timing report:", str(data_list))
    
    fields = QgsFields() #there's a BUG in QGIS here (?), normally : fields = ....
    fields.append(QgsField("Source", QVariant.String ))
    fields.append(QgsField("Target", QVariant.String))
## fields.append(QgsField("Source_lbl", QVariant.String, 'string',50))
## fields.append(QgsField("Target_lbl", QVariant.String, 'string',50))
    fields.append(QgsField("Visible", QVariant.String, 'string',5))
    fields.append(QgsField("TargetSize", QVariant.Double, 'double',10,3))
    fields.append(QgsField("Distance", QVariant.Double, 'double',10,2))

    writer = QgsVectorFileWriter( file_name + ".shp", "CP1250", fields,
                                  QGis.WKBLineString, coordinate_ref_system) #, "ESRI Shapefile"
                                            #CP... = encoding
    if writer.hasError() != QgsVectorFileWriter.NoError:
        QMessageBox.information(None, "ERROR!", "Cannot write intervisibilty file (?)")
        return 0
    
    for r in data_list:
        # create a new feature
        feat = QgsFeature()
        if use_pix_coords:
            half_pix= pix/2 #global variable pix
            l_start=QgsPoint(raster_x_min  + r[1]*pix + half_pix, raster_y_max - r[2]*pix - half_pix )
            l_end = QgsPoint(raster_x_min  + r[4]*pix + half_pix, raster_y_max - r[5]*pix - half_pix)
        else:
            l_start=QgsPoint(r[1],r[2]);  l_end = QgsPoint(r[4],r[5])
         
        feat.setGeometry(QgsGeometry.fromPolyline([l_start, l_end]))
        # do not cast ID to string: unicode problem -- angle * distance in pixels -- distance * pixel_size
        #feat.setAttributes([ str(r[0]), str(r[3]), bool(r[6]), float(r[7] * r[8]), ])
        feat.setFields(fields)
        feat['Source'] = r[0]
        feat['Target'] = r[3]
        feat['Visible'] = 'True' if r[6] else 'False'
        feat['TargetSize'] = float(r[7])
        feat['Distance'] = float(r[8])
        
        writer.addFeature(feat)
        del feat

    del writer
    layer = None
    return file_name + ".shp"


