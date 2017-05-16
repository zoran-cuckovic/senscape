# -*- coding: utf-8 -*-

#from qgis.core import Qgis, QgsUnitTypes
from PyQt4.QtGui import QMessageBox
from qgis.core import *
from os import path

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
from processing.core.parameters import (ParameterVector,
                                        ParameterRaster,
                                        ParameterNumber,
                                        ParameterBoolean,
                                        ParameterSelection,
                                        ParameterTableField,
                                        ParameterFile)
from processing.core.outputs import Output, OutputRaster
from processing.tools import dataobjects
from processing.core.ProcessingLog import ProcessingLog
from processing.gui import MessageBarProgress

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
        
       
        self.addOutput(OutputRaster(
            self.OUTPUT,
            self.tr('Output file')))
    
    def help(self):
        return False, 'http://zoran-cuckovic.github.io/senscape/help/raster'
        

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



        #getTempFilenameInTempFolder(
            #self.name + '.' + self.getDefaultFileExtension(alg)

            
        # output_dir = self.getOutputValue(self.OUTPUT_DIR)

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
                         
        points = pts.Points(observers)

        miss = points.test_fields(["observ_hgt", "radius"])

        if miss:
            QMessageBox.information(None, "ERROR!",
                "Missing fields! \n" + "\n".join(miss))
            return


        points.take(dem.extent, dem.pix)

        if points.count == 0:
            QMessageBox.information(None, "ERROR!",
                "No viewpoints in the chosen area!")
            return
        
        #special case: singe file
        elif points.count == 1: operator = -1
        
        # second test, after checking the number of points       
        elif operator == 0:
            miss = points.test_fields(["file"])
            if miss:
                QMessageBox.information(None, "ERROR!",
                "Missing fields : file !")
                return

  
            
            
        dem.set_buffer_mode(operator)

        prog = MessageBarProgress.MessageBarProgress()

	#will assign the result to dem class (?)	
        files, report =  ws.Viewshed (points, dem, 
                                analysis_type,
                              Target_points=None,
                              curvature=useEarthCurvature,
                                refraction=refraction,
                                algorithm = precision,
                                progress = prog )

        prog.close()

        
        if operator==0:
            #this is a hack: to clear outputs
            #(there is a function "removeOutputFromName" but this is simpler ..) 
            self.outputs=[]
            for f in files:
                #self.outputs = list of instances !
                o = OutputRaster( "dd")
                o.setValue(f)
                # there is also .name property, but it's not working (?!)
                #o.description = path.basename(f)[:-4]

                self.addOutput(o)

                #name = os.path.basename(new_file)
                #ProcessingResults.addResult(name, new_file)

                # this is set in proseccing options !!
##        else:
##             self.outputs[0].description = path.basename(output_path)[:-4]
            

        txt = "Analysed points \n ID : visible pixels : total area" 
        
        for l in report:
            txt = txt + "\n" + ' , '.join(str(x) for x in l)
            
        ProcessingLog.addToLog ('INFO', txt) 

        
        """
        Is this a good practice? dem object has been modified in ws.Viewshed routine.
        But there is no sense in putting viewshed algorithm in Raster class,
        it is also used for vectors (intervisibility). And it's a lot of code .. 

        """
