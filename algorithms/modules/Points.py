# -*- coding: utf-8 -*-

from PyQt4.QtCore import *
from qgis.core import *

from os import path

import gdal #ogr

import numpy as np

import Raster as rst


"""
Points class is creating a clean shapefile with analysis parameteres in
the associated table. It is also taking care of geometries (filtering
points outside raster, searching for neighbours in a specified radius etc.)
The idea is to move all the mess of handling vector input and output
in a single class.
"""

class Points:
    
    def __init__(self, shapefile):
    #layer = qgis layer,  bounding_box = QgsRectangle,   field_id= string)
        self.pt={}

        self.layer = QgsVectorLayer(shapefile, 'o', 'ogr')

        self.crs = self.layer.crs()

        #make a test !
        self.missing = []

        fields = ["ID", "observ_hgt",  "radius"]
        provider = self.layer.dataProvider()

     
        for f in fields:
            if provider.fieldNameIndex(f)==-1:
                self.missing.append(f)

        # optional missing
        
        
        self.count = 0 # only take routine can determine the number of used points
        

    """
    Take care of parameters. Leave geometries for other routine:
    same points can be used with different rasters. 

    """

    # multiple parameters (numeric + field name) enable to set a default
    # in case there is a problem, eg z = float(None)
    def clean_parameters(self, z_obs, radius ,
                        z_targ=0, 
                        field_ID = None,
                        field_zobs = None,
                        field_ztarg=None,
                        field_radius=None,
                        folder = None):
       

        errors=[]
       
        
        #feat = QgsFeature()
        
        #feature_iterator= layer.getFeatures() 
        
        #while feat_iterator.nextFeature(feat):

        for feat in self.layer.getFeatures():

            
            
            geom = feat.geometry()
            t = geom.asPoint()

            x_geog, y_geog= t

            z,zt,r = z_obs, z_targ, radius
         

            if field_ID: id1 = feat[field_ID]
            else : id1 = feat.id()

            # test for duplicates
            if id1 not in self.pt:
                
                #addition for possible field values.
                #override with fixed parameters in case of problem 
                if field_zobs :
                    try : z = float(feat[field_zobs])
                    except: z=z_obs

                if field_radius:
                    try : r = float(feat[field_radius]) 
                    except: r=radius
                    
                # obligatory prarameters        
                self.pt[id1]={"z":z ,  "radius" : r,
                              "x_geog":x_geog, "y_geog" : y_geog }

                # optional 

                if field_ztarg:
                    try : self.pt[id1]["z_targ"] = float(feat[field_target])
                    except: self.pt[id1]["z_targ"]=z_targ

                if folder:
                    self.pt[id1]["path"] = path.join(folder, str(id1) + ".tif")

                    
     
            else: errors.append(id1)

        
        return errors if errors else 0
         
       # self.max_radius = max(x, key=lambda i: x[i])

    """
        
    Find the highest point in a perimeter around each observer point.
    Note that it always moves the point to the center of the highest pixel

    This is duplicating functions from take_points (selection inside a frame)!
    to reorganise !! (should be independent function, could be useful elsewhere...)
       
    """
    def move_top(self,  raster_path, search_radius):

  
        r = rst.Raster(raster_path)
        pix = r.pix; half_pix = pix/2

        raster_x_min, raster_y_min, raster_x_max,raster_y_max = r.extent
                
        radius_pix = int(search_radius/pix)

        win_size = radius_pix * 2 + 1

        r.set_master_window(radius_pix)
        
        
    #raster_y_min = raster_y_max - raster_y_size * pix
    #raster_x_max = raster_x_min + raster_x_size * pix

        for key in self.pt:

                     
            x, y = self.pt[key]["x_geog"], self.pt[key]["y_geog"]

            #how to use Qgs functions?
            # ext = QgsRectangle(*extents)
            #pt = QgsPoint(x,y)
            # pt in ext ???

          
            if not raster_x_min < x < raster_x_max \
            or not raster_y_min < y < raster_y_max: continue


            # make a function for pixel coords !!!
            x_pix= int((x - raster_x_min) / pix)
            y_pix = int((raster_y_max - y) / pix) #reversed !

            r.open_window(x_pix, y_pix, radius_pix, initial_value=r.min)

            #chunks are padded to be square (for viewsheds)
            # x, y is always in the centre
            
            iy, ix=np.unravel_index( np.argmax(r.window),
                                     r.window.shape )
        
            # unravel is giving offsets inside the window,
            # we need to place it inside the entire raster
            x_off , y_off, win_x, win_y = r.gdal_slice
            # when the point is close to border, take into account window overlap!
            x_off -= win_size - win_x; y_off -= win_size - win_y

            self.pt[key]["x_geog"] = (ix + x_off) * pix + raster_x_min + half_pix
            self.pt[key]["y_geog"] = raster_y_max - (iy + y_off) * pix  - half_pix

#           
##            if iy != radius_pix or ix != radius_pix:
##        
##                self.pt[key]["x_geog"] += (ix - x_pix)  * pix
##                self.pt[key]["y_geog"] += (iy - y_pix)  * pix
        

    """
                
            # we cannot know the position of the observer! if it is not in the center ...
            z_top = None
            
            for j in xrange(0, y_size): 
                for i in xrange(0, x_size):
                    try: k = dt [j, i] # it may be an empty cell or whatever...
                    except: continue
                    
                    if k > z_top: x_top,y_top,z_top = i,j,k

            if x_off1: x_top = pt_x + (x_top - search_top)
            if y_off1: y_top = pt_y + (y_top - search_top)



            #todo                 
            x_geog += (x2 - x)  * pix
            y_geog += (y2 - y) * pix
    """

    
    """ much faster with numpy ...

    if not 0 <= x < raster_x_size or not 0 <= y < raster_y_size : continue
       
       #cropping from the front
    if x <= radius_pix:   x_offset =0
       #cropping from the back
    else:      x_offset = x - radius_pix         
      
    if y <= radius_pix:  y_offset =0
    else:   y_offset = y - radius_pix

    x_offset2 = min(x + radius_pix +1, raster_x_size) #could be enormus radius, so check both ends always
    y_offset2 = min(y + radius_pix + 1, raster_y_size )
    
    window_size_y = y_offset2 - y_offset
    window_size_x = x_offset2 - x_offset

    mx = r.ReadAsArray(x_offset, y_offset, window_size_x, window_size_y).astype(float)
    m = np.argmax(mx)

    iy, ix=np.unravel_index(m, mx.shape)
    
    #0.5 is to move to the center of corresp. pixel
    x2_g = ( ix +0.5 + x_offset) * pix  + raster_x_min 
    
    y2_g = raster_y_max - (y_offset + iy + 0.5) * pix  

    g= QgsGeometry.fromPoint(QgsPoint( x2_g, y2_g))
    
    inputLayer.dataProvider().changeGeometryValues({ pt_id : g })

    """
    """
    Assign targets for each observer point. Targets are a class instance
    Used after take which selects good points
    """
    def network (self,targets):

        self.edges={}

        for pt1 in self.pt: 

                      
            x = self.pt[pt1]["x_pix"]
            y = self.pt[pt1]["y_pix"]
            r = self.pt[pt1]["radius"] #it's pixelised after take !!

            radius_pix= int(r)
            max_x = x + radius_pix; min_x = x - radius_pix
            max_y = y + radius_pix; min_y = y - radius_pix
            #does not need cropping if target points match raster extent


            # cheap distance check 
            #   if mx_dist [y2_local,x2_local] > radius: continue
            #local coords :in intervisibilty

                
##                if z_target_field: #this is a clumsy addition so that each point might have it's own height
##                    try: tg_offset = float(feat2[z_target_field])
##                    except: pass
            
            for pt2 in targets.pt:
                x2 = targets.pt[pt2]["x_pix"]
                y2 = targets.pt[pt2]["y_pix"]

                #skipping 1
                if x2==x and y2==y : continue
                
                if min_x <= x2 <= max_x and min_y <= y2 <= max_y:
                    if  (x-x2)**2 + (y-y2)**2 <= r**2:
                          # this is inefficient for looping
                          # need to open a window for each edge...
##                        self.edges[id1,id2]={}
                        try: self.pt[pt1]["targets"].append(pt2)
                        except: self.pt[pt1]["targets"]=[pt2]


    """
    Returns a dict of points, prepared for visibilty analysis
    All values are expressed in pixel offsets.

    This must work on a freshly loaded shapefile,
    where points.missing is [] !!
    [_or make a new function to test shapefiles_]
   
    """
    def take (self, extent, pix_size, spatial_index=None):
        
        self.max_radius = 0

        x_min, y_max = extent[0], extent [3]

        bounding_box = QgsRectangle(*extent) #* unpacks an argument list     

        if not spatial_index: #for intersect, not very helpful ...?
            s_index = QgsSpatialIndex()
            for f in self.layer.getFeatures():  s_index.insertFeature(f)
            
        else: s_index = spatial_index

        feature_ids = s_index.intersects(bounding_box)

        for fid in feature_ids:
            feat = self.layer.getFeatures(
                   QgsFeatureRequest().setFilterFid(fid)).next()

            geom = feat.geometry()
            t = geom.asPoint()

            x_geog, y_geog= t
            
            x = int((x_geog - x_min) / pix_size) # not float !
            y = int((y_max - y_geog) / pix_size) #reversed !

            r = feat["radius"] / pix_size

            if r > self.max_radius : self.max_radius = r
            
            try: tg = feat["target_hgt"]
            except : tg = 0

            try: f = feat["file"]
            except: f = None
            
            self.pt[ feat["ID"] ]={"z" : feat["observ_hgt"] ,
                                    "z_targ": tg,
                                    "radius" : r,
                                    "x_pix" : x, "y_pix":y,
                                    "file" : f}       
                                    #"x_geog" :x_geog, "y_geog": y_geog,
        
        self.count = len(feature_ids)
        


                        
    def write_points (self, file_name, coordinate_ref_system, use_pix_coords=False):

        from processing.tools.vector import VectorWriter


        # test the existence of values in dict 
        miss_tg = "target_hgt" not in self.pt.values()[0]
        miss_file= "path" not in self.pt.values()[0]
        
        
        fields = QgsFields()
        fields.append ( QgsField("ID", QVariant.String, 'string',50))
        fields.append (QgsField("observ_hgt", QVariant.Double,'double', 5,4 ))
        
        fields.append ( QgsField("radius", QVariant.Double, 'double',10,3))
        if not miss_tg:
            fields.append (QgsField("target_hgt", QVariant.Double,'double', 5,4 ))
        if not miss_file:
            #windows has 260 character limit for filenames...
            fields.append ( QgsField("file", QVariant.String, 'string',260))

        

       # writer = QgsVectorFileWriter( file_name + ".shp", "CP1250", fields,
       #                               QGis.WKBPoint, coordinate_ref_system, "ESRI Shapefile")

      
        writer = VectorWriter(file_name, None, fields, QGis.WKBPoint,
                              coordinate_ref_system)
        
                                                #CP... = encoding
##        if writer.hasError() != QgsVectorFileWriter.NoError:
##            QMessageBox.information(None, "ERROR!", "Cannot write points file (?)")
##            return 0
        
        for r in self.pt :

           
            # create a new feature
            feat = QgsFeature()
            
            feat.setGeometry(QgsGeometry.fromPoint(
                QgsPoint(float(self.pt[r]["x_geog"]),
                         float(self.pt[r]["y_geog"] )) ))

      
            feat.setFields(fields)

            feat['ID'] = r
            feat ['observ_hgt']=self.pt[r]["z"]
            feat ['radius']=self.pt[r]["radius"]
            if not miss_tg:
                feat ['target_hgt']=self.pt[r]["z_targ"]
            if not miss_file:
                feat ['file']=self.pt[r]["path"]
            writer.addFeature(feat)

        del writer

        
       
        #return file_name + ".shp"



    def write_network (self, file_name, edges, target_class,
                   coordinate_ref_system=None , use_pix_coords=False):

        from processing.tools.vector import VectorWriter

        #QMessageBox.information(None, "Timing report:", str(data_list))
        fields = QgsFields()
        fields.append ( QgsField("Source", QVariant.String, 'string',50))
        fields.append ( QgsField("Target", QVariant.String, 'string',50))
       # fields.append (QgsField("Visible", QVariant.Boolean ))
        fields.append(QgsField("Visible", QVariant.String, 'string',5))
        fields.append (QgsField("observ_hgt", QVariant.Double,'double', 5,4 ))
        fields.append (QgsField("target_hgt", QVariant.Double,'double', 5,4 ))
        

        

       # writer = QgsVectorFileWriter( file_name + ".shp", "CP1250", fields,
       #                               QGis.WKBPoint, coordinate_ref_system, "ESRI Shapefile")
        
        tg= target_class

        crs = coordinate_ref_system if coordinate_ref_system else self.crs
        
        writer = VectorWriter(file_name, None, fields, QGis.WKBLineString,crs)
        
                                                #CP... = encoding
##        if writer.hasError() != QgsVectorFileWriter.NoError:
##            QMessageBox.information(None, "ERROR!", "Cannot write points file (?)")
##            return 0
        
        for r in self.pt :

            if "targets" not in self.pt[r]: continue
           
            # create a new feature
            feat = QgsFeature()

            p1 = QgsPoint(float(self.pt[r]["x_geog"]),
                         float(self.pt[r]["y_geog"] ))

            for t in self.pt[r]["targets"]:

                p2 = QgsPoint(float(tg.pt[t]["x_geog"]),
                         float(tg.pt[t]["y_geog"] ))
            
    
                feat.setGeometry(QgsGeometry.fromPolyline([p1, p2]))
          
                feat.setFields(fields)

                feat['Source'] = r
                feat['Target'] = t
                feat['Visible'] = 'True' if edges[r,t] >= 0 else 'False'
                feat ['observ_hgt']=self.pt[r]["z"]
                feat ['target_hgt']=float(edges[r,t]) #why doesn't it accept Python numbers ??        
                
                writer.addFeature(feat)

        del writer
        

    
