# -*- coding: utf-8 -*-

#from qgis.core import Qgis, QgsUnitTypes
from PyQt4.QtGui import QMessageBox

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
from processing.core.parameters import (ParameterVector,
                                        ParameterRaster,
                                        ParameterNumber,
                                        ParameterBoolean,
                                        ParameterSelection,
                                        ParameterTableField)
from processing.core.outputs import OutputRaster, OutputDirectory
from processing.tools import dataobjects
from processing.core.ProcessingLog import ProcessingLog

from .modules import doViewshed as ws
from .modules import Points as pts
from .modules import Raster as rst


class Viewshed(GeoAlgorithm):

    DEM = 'DEM'
    OBSERVER_POINTS = 'OBSERVER_POINTS'
    
    USE_CURVATURE = 'USE_CURVATURE'
    REFRACTION = 'REFRACTION'
    PRECISION = 'PRECISION'
    ANALYSIS_TYPE = 'ANALYSIS_TYPE'
    OPERATOR = 'OPERATOR'
    OUTPUT = 'OUTPUT'
    OUTPUT_DIR = 'OUTPUT_DIR'

    PRECISIONS = ['Coarse','Normal', 'Fine']
    
    TYPES = ['Binary viewshed', 'Depth below horizon',
             'Horizon', 'Horizon - intermediate']

    OPERATORS = ['Individual files', 'Sum', "Maximum", "Minimum"]

    def defineCharacteristics(self):
        self.name = 'Viewshed'
        self.group = 'Visibility analysis'


        self.addParameter(ParameterVector(
            self.OBSERVER_POINTS,
            self.tr('Observer location(s)'),
            0))  #dataobjects.TYPE_VECTOR_POINT

        self.addParameter(ParameterRaster(
            self.DEM,
            self.tr('Digital elevation model')))

        self.addParameter(ParameterSelection(
            self.ANALYSIS_TYPE,
            self.tr('Analysis type'),
            self.TYPES,
            0))

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
            self.OPERATOR,
            self.tr('Combining multiple outputs'),
            self.OPERATORS,
            0))
        # would be cool if this widget get locked/unlocked according to chosen operator
        self.addOutput(OutputDirectory(
            self.OUTPUT_DIR,
            self.tr('Output directory (for individual output)')))
        # OutputDirectory( ...
        self.addOutput(OutputRaster(
            self.OUTPUT,
            self.tr('Output file')))

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

        
        useEarthCurvature = self.getParameterValue(self.USE_CURVATURE)
        refraction = self.getParameterValue(self.REFRACTION)
        precision = self.getParameterValue(self.PRECISION)
        analysis_type = self.getParameterValue(self.ANALYSIS_TYPE)
        operator =  self.getParameterValue(self.OPERATOR)
        

        output_path = self.getOutputValue(self.OUTPUT)
        output_dir = self.getOutputValue(self.OUTPUT_DIR)
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



      
                    
        dem = rst.Raster(raster, output=output_path)
        # ADD SOME TESTS: points out of raster [OK],
        # raster rotated [projections ??], rectnagular pixels [OK?]
        """
        In 'Raster' class : self.crs = crs if crs else gdal_raster.GetProjection()
        
        Perhaps there could be problems when registered crs
        is not matching the one chosen in QGIS ??
        
        ... to override one can do dem.crs = ....
        """       
                         
        points = pts.Points(observers)#
        #normally  .missing is an (empty) list
        if  points.missing :
            
            QMessageBox.information(None, "ERROR!",
            "Missing fields! \n" + "\n".join(points.missing))
            return

        points.take(dem.extent, dem.pix)


        if points.count == 0: pass # RAISE ERROR!

        #special case: singe file
        elif points.count == 1: operator = -1

        else:
            dem.set_buffer(operator,
                           fill_nan = True if operator > 1
                                            else False) #default is 0

        # to avoid passing this argument to doViewshed, register as property
        # the idea is to clean doViewshed from handling raster input/output (??)
        # (normally output path is passed upon creation of the class)
        if operator == 0:
            dem.output = output_dir

	#will assign the result to dem class (?)	
        report =  ws.Viewshed (points, dem, 
                                analysis_type,
                              Target_points=None,
                              curvature=useEarthCurvature,
                                refraction=refraction,
                                algorithm = 1)

        txt = "Analysed points \n ID : visible pixels" 
        
        for l in report:
            txt = txt + "\n" + ' , '.join(str(x) for x in l)
            
        ProcessingLog.addToLog ('INFO', txt) 

        

        """
        Is this a good practice? dem object has been modified in ws.Viewshed routine.
        But there is no sense in putting viewshed algorithm in Raster class,
        it is also used for vectors (intervisibility). And it's a lot of code .. 

        """
