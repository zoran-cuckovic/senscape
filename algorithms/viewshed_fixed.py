# -*- coding: utf-8 -*-

#from qgis.core import Qgis, QgsUnitTypes

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
from processing.core.parameters import (ParameterVector,
                                        ParameterRaster,
                                        ParameterNumber,
                                        ParameterBoolean,
                                        ParameterSelection,
                                        ParameterTableField)
from processing.core.outputs import OutputRaster
from processing.tools import dataobjects

from .modules import doViewshed as ws
from .modules import Points as pts
from .modules import Raster as rst


class ViewshedFixed(GeoAlgorithm):

    DEM = 'DEM'
    OBSERVER_POINTS = 'OBSERVER_POINTS'
    OBSERVER_ID = 'OBSERVER_ID'
    OBSERVER_HEIGHT = 'OBSERVER_HEIGHT'
    TARGET_HEIGHT = 'TARGET_HEIGHT'
    SEARCH_RADIUS = 'SEARCH_RADIUS'
    USE_CURVATURE = 'USE_CURVATURE'
    REFRACTION = 'REFRACTION'
    PRECISION = 'PRECISION'
    ANALYSIS_TYPE = 'ANALYSIS_TYPE'
    CUMULATIVE = 'CUMULATIVE'
    OUTPUT = 'OUTPUT'

    PRECISIONS = [ 'Normal']  # 'Coarse',   , 'Fine']
    TYPES = ['Binary viewshed'] #, 'Invisibility depth', 'Horizon', 'Horizon full']

    def defineCharacteristics(self):
        self.name = 'Viewshed (fixed)'
        self.group = 'Visibility analysis'

        self.addParameter(ParameterRaster(
            self.DEM,
            self.tr('Digital elevation model')))
        self.addParameter(ParameterVector(
            self.OBSERVER_POINTS,
            self.tr('Observer location(s)'),
            0))  #dataobjects.TYPE_VECTOR_POINT
##        self.addParameter(ParameterTableField(
##            self.OBSERVER_ID,
##            self.tr('Observer ids (leave unchanged to use feature ids)'),
##            self.OBSERVER_POINTS,
##            optional=True))
        self.addParameter(ParameterNumber(
            self.OBSERVER_HEIGHT,
            self.tr('Observer height, meters'),
            0.0, 999.999, 1.6))
        self.addParameter(ParameterNumber(
            self.TARGET_HEIGHT,
            self.tr('Target height, meters'),
            0.0, 999.999, 0.0))
        self.addParameter(ParameterNumber(
            self.SEARCH_RADIUS,
            self.tr('Search radius, meters'),
            0.0, 99999999.99999, 5000))
        self.addParameter(ParameterBoolean(
            self.USE_CURVATURE,
            self.tr('Take in account Earth curvature'),
            False))
        self.addParameter(ParameterNumber(
            self.REFRACTION,
            self.tr('Atmoshpheric refraction'),
            0.0, 1, 0.13))
        self.addParameter(ParameterSelection(
            self.PRECISION,
            self.tr('Algorithm precision'),
            self.PRECISIONS,
            1))
        self.addParameter(ParameterSelection(
            self.ANALYSIS_TYPE,
            self.tr('Analysis type'),
            self.TYPES,
            0))
        # this is by default for fixed versions
##        self.addParameter(ParameterBoolean(
##            self.CUMULATIVE,
##            self.tr('Generate cumulative output'),
##            False))

        # OutputDirectory( ...
        self.addOutput(OutputRaster(
            self.OUTPUT,
            self.tr('Output')))

    def processAlgorithm(self, feedback):

        # objects are not used : the idea is to be able to use modules outside processing
        # it's simpler to pass the path (?)
##        raster_object = dataobjects.getObjectFromUri(
##            self.getParameterValue(self.DEM))
##
##        observer_object = dataobjects.getObjectFromUri(
##            self.getParameterValue(self.OBSERVER_POINTS))

        raster = self.getParameterValue(self.DEM)
        observers = self.getParameterValue(self.OBSERVER_POINTS)
#        observerIdField = self.getParameterValue(self.OBSERVER_ID)
        observer_height = self.getParameterValue(self.OBSERVER_HEIGHT)
        target_height = self.getParameterValue(self.TARGET_HEIGHT)
        search_radius = self.getParameterValue(self.SEARCH_RADIUS)
        useEarthCurvature = self.getParameterValue(self.USE_CURVATURE)
        refraction = self.getParameterValue(self.REFRACTION)
        precision = self.getParameterValue(self.PRECISION)
        analysis_type = self.getParameterValue(self.ANALYSIS_TYPE)
        cumulative = self.getParameterValue(self.CUMULATIVE)

        output_path = self.getOutputValue(self.OUTPUT)

        # convert meters to layer distance units
        # [this can be confusing when the module is used in a script,
        #  and it's 3.0 function ]
        #coef = QgsUnitTypes.fromUnitToUnitFactor(Qgis.DistanceMeters, dem.crs().mapUnits())
		
        #searchRadius = searchRadius * coef


##        visibility.viewshed(raster,
##                            observer,
##                            observerIdField,
##                            observerHeight,
##                            targetHeight,
##                            searchRadius,
##                            useEarthCurvature,
##                            refraction,
##                            analysisType,
##                            precision,
##                            outputPath)



      
                        
        dem = rst.Raster(raster)
        # crs=  layer.crs().toWkt() #I do not have layer object ??
        """
        In 'Raster' class : self.crs = crs if crs else gdal_raster.GetProjection()
        
        Perhaps there could be problems when registered crs
        is not matching the one chosen in QGIS ??
        
        ... to override one can do Dem_raster.crs = ....
        """

       
                         
        points = pts.Points(observers, dem.extent, dem.pix,
                            observer_height, search_radius,
                            target_height) # and all other stuff ....

          # to add up results
        if points.count >1:
            #cannot determine size, only write mode : to implement in the future...
            dem.set_buffer_mode ( 1)

        
        
	#will assign the result to dem class (?)

        	
        success =  ws.Viewshed (points, dem, 
                                analysis_type,
                              Target_points=None,
                              curvature=useEarthCurvature,
                                refraction=refraction,
                                algorithm = 1)

        

        dem.write_result(output_path)

        """
        Is this a good practice? dem object has been modified in ws.Viewshed routine.
        But there is no sense in putting viewshed algorithm in Raster class,
        it is also used for vectors (intervisibility). And it's a lot of code .. 

        """
