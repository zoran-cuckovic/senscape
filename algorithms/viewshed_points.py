# -*- coding: utf-8 -*-

"""
/***************************************************************************
 TestProcessing
                                 A QGIS plugin
 Some descr
                              -------------------
        begin                : 2017-03-10
        copyright            : (C) 2017 by some
        email                : some
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'some'
__date__ = '2017-03-10'
__copyright__ = '(C) 2017 by some'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from PyQt4.QtCore import QSettings
from qgis.core import QgsVectorFileWriter

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.parameters import (ParameterRaster,
                                        ParameterVector,
                                        ParameterBoolean,
                                        ParameterString,
                                        ParameterSelection,
					ParameterTableField,
                                        ParameterNumber,
                                        ParameterRaster)	
from processing.core.outputs import  OutputVector
from processing.tools import dataobjects, vector


from PyQt4.QtGui import * #this is for message box : should use processing log ?


class ViewshedPoints(GeoAlgorithm):
    """This is an example algorithm that takes a vector layer and
    creates a new one just with just those features of the input
    layer that are selected.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the GeoAlgorithm class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OBSERVER_POINTS = 'OBSERVER_POINTS'
    INPUT_DEM = 'INPUT_RASTER'

    OUTPUT_VECTOR = 'OUTPUT_VECTOR'

    OBSERVER_ID = 'OBSERVER_ID'

    RADIUS = 'RADIUS'
    RADIUS_FIELD = 'RADIUS_FIELD'
    
    OBS_HEIGHT = 'OBS_HEIGHT'
    OBS_HEIGHT_FIELD = 'OBS_HEIGHT_FIELD'

    TARGET_HEIGHT = 'TARGET_HEIGHT'
    TARGET_HEIGHT_FIELD = 'TARGET_HEIGHT_FIELD'

    MOVE_TOP = 'MOVE_TOP'
    
    def defineCharacteristics(self):
        """Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # The name that the user will see in the toolbox
        self.name = 'Create viewpoints'

        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Visibility analysis'

        # We add the input vector layer. It can have any kind of geometry
        # It is a mandatory (not optional) one, hence the False argument    


        self.addParameter(ParameterVector(
            self.OBSERVER_POINTS,
            self.tr('Observer location(s)'),
            0)) #dataobjects.TYPE_VECTOR_POINT
       
        self.addParameter(ParameterTableField(
            self.OBSERVER_ID,
            self.tr('Observer ids (leave unchanged to use feature ids)'),
            self.OBSERVER_POINTS,
            optional=True))

        self.addParameter(ParameterNumber(
            self.RADIUS,
            self.tr("Radius of analysis, meters"),
            0.0, 99999999.999, 5000))
  
        self.addParameter(ParameterTableField(
            self.RADIUS_FIELD,
            self.tr('Field value for analysis radius'),
            self.OBSERVER_POINTS,
            optional=True))

        self.addParameter(ParameterNumber(
            self.OBS_HEIGHT,
            self.tr('Observer height, meters'),
            0.0, 999.999, 1.6))
        
        self.addParameter(ParameterTableField(
            self.OBS_HEIGHT_FIELD,
            self.tr('Field value for observer height'),
            self.OBSERVER_POINTS,
            optional=True))

        self.addParameter(ParameterNumber(
            self.TARGET_HEIGHT,
            self.tr('Target height, meters'),
            0.0, 999.999, 0.0))
        
        self.addParameter(ParameterTableField(
            self.TARGET_HEIGHT_FIELD,
            self.tr('Target height, meters'),
            self.OBSERVER_POINTS,
            optional=True))

        self.addParameter(ParameterNumber(
            self.MOVE_TOP,
            self.tr('Find highest elevation, radius in meters'),
            0.0, 9999.99, 0.0))
##
        self.addParameter(ParameterRaster(self.INPUT_DEM,
            self.tr('Elevation model for moving points'),  False))


        self.addOutput(OutputVector(self.OUTPUT_VECTOR,
                                    self.tr('Output viewshed points')))
# optional param
#        self.addParameter(ParameterString(self.RADIUS,
#                                          self.tr("SOME MOCK PARAMETER"),
#                                                    '', optional=True)) 
        # We add a vector layer as output
      ##    NEW ?
##                self.addParameter(ParameterSelection(self.RTYPE,
##                                             self.tr('Output raster type'),
##                self.TYPE, 5))

    def processAlgorithm(self, progress):
        
	
	# is it better to import here or on top ?
        import gdal
        from qgis.core import QgsRectangle
        import numpy as np
        
        
        from .modules import Points as pts
        from .modules import Raster as rst

        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        #Raster_path = self.getParameterValue(self.INPUT_DEM)
        Points_path = self.getParameterValue(self.OBSERVER_POINTS)
        
        Output = self.getOutputValue(self.OUTPUT_VECTOR)

        observer_id = self.getParameterValue(self.OBSERVER_ID)
        
        observer_height = self.getParameterValue(self.OBS_HEIGHT)
        observer_height_field =  self.getParameterValue(self.OBS_HEIGHT_FIELD)
        
        radius = self.getParameterValue(self.RADIUS)
        radius_field = self.getParameterValue(self.RADIUS_FIELD)
        
        target= self.getParameterValue(self.TARGET_HEIGHT)
        target_field= self.getParameterValue(self.TARGET_HEIGHT_FIELD)

        move = self.getParameterValue(self.MOVE_TOP)

        
        
                        
        
                        # crs=  layer.crs().toWkt() #I do not have layer object ??
                        # write_mode = 'cumulative', if cumulative
        
        
        points = pts.Points(Points_path) # and all other stuff ....

        success = points.clean_parameters( observer_height, radius,
                           z_targ = target ,
                           field_ID = observer_id,
                           field_zobs = observer_height_field,
                           field_ztarg=target_field,
                           field_radius=radius_field)
        
        if success != 0 :
            print success
            QMessageBox.information(None, "Duplicate IDs!", str(success))
            return


        if move:
            

            points.move_top(self.getParameterValue(self.INPUT_DEM), move)
        
        points.write_points (Output, points.crs)
        
        """ 
        Perhaps there could be problems when registered crs
        is not matching the one chosen in QGIS ??
        In 'Raster' class : self.crs = crs if crs else gdal_raster.GetProjection()
        so to override on can do Dem_raster.crs = ....
        """
       # dem.write_result(Output)
        

        """
        Is this a good practice? dem object has been modified in ws.Viewshed routine.
        But there is no sense in putting viewshed algorithm in Raster class,
        it is also used for vectors (intervisibility). And it's a lot of code .. 

        """

        # write output
        #driver = gdal.GetDriverByName( 'GTiff' )
          
        #c = driver.CreateCopy(Output, r , 0) #WHAT IS 0 ? : "strict or not" default =1

        #c.GetRasterBand(1).WriteArray(matrix_final)

        #c = None
