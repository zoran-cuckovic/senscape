# -*- coding: utf-8 -*-

#from qgis.core import Qgis, QgsUnitTypes
from PyQt4.QtGui import QMessageBox
from os import path

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
from processing.core.parameters import (ParameterVector,
                                        ParameterRaster,
                                        ParameterNumber,
                                        ParameterBoolean,
                                        ParameterSelection,
                                        ParameterTableField)
from processing.core.outputs import OutputVector
from processing.tools import dataobjects


from .modules import doViewshed as ws
from .modules import Points as pts
from .modules import Raster as rst


class Intervisibility(GeoAlgorithm):

    DEM = 'DEM'
    
    OBSERVER_POINTS = 'OBSERVER_POINTS'
    TARGET_POINTS = 'TARGET_POINTS'
    
    USE_CURVATURE = 'USE_CURVATURE'
    REFRACTION = 'REFRACTION'

    OUTPUT_VECTOR = 'OUTPUT_VECTOR'


    PRECISIONS = ['Coarse', 'Normal']
    PRECISION = 'PRECISION'
    

    def defineCharacteristics(self):
        self.name = 'Intervisibility'
        self.group = 'Visibility analysis'

        self.addParameter(ParameterRaster(
            self.DEM,
            self.tr('Digital elevation model')))

        self.addParameter(ParameterVector(
            self.OBSERVER_POINTS,
            self.tr('Observer location(s)'),
            0)) #dataobjects.TYPE_VECTOR_POINT

        self.addParameter(ParameterVector(
            self.TARGET_POINTS,
            self.tr('Target location(s)'),
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

        self.addOutput(OutputVector(
            self.OUTPUT_VECTOR,
            self.tr('Output')))

    def processAlgorithm(self, feedback):
##        dem = dataobjects.getObjectFromUri(
##            self.getParameterValue(self.DEM))
##        observer = dataobjects.getObjectFromUri(
##            self.getParameterValue(self.OBSERVER_POINTS))

        raster = self.getParameterValue(self.DEM)
        observers = self.getParameterValue(self.OBSERVER_POINTS)
        targets = self.getParameterValue(self.TARGET_POINTS)
        
        useEarthCurvature = self.getParameterValue(self.USE_CURVATURE)
        refraction = self.getParameterValue(self.REFRACTION)
        precision = self.getParameterValue(self.PRECISION)
   
        output_path = self.getOutputValue(self.OUTPUT_VECTOR)

        # convert meters to layer distance units
##        coef = QgsUnitTypes.fromUnitToUnitFactor(Qgis.DistanceMeters, dem.crs().mapUnits())
##
##
##        searchRadius = searchRadius * coef
##
##
##        visibility.viewshed(dem,
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

       
        o= pts.Points(observers)
        
        t= pts.Points(targets)

        required=["observ_hgt", "radius"]

        miss1 = o.test_fields (required)
        miss2 = t.test_fields (required)

        if miss1 or miss2:
            QMessageBox.information(None, "Missing fields!",
                "Missing in observer points: " + str(miss1) +
                "\n Missing in target points: " + str(miss2))
            return

        o.take(dem.extent, dem.pix)
        t.take(dem.extent, dem.pix)

        o.network(t) #do this after .take which takes points within raster extents
        	
        relations =  ws.intervisibility2(o,t, dem,                                   
                                    curvature=useEarthCurvature,
                                    refraction=refraction,
                                    interpolate = precision)

        o.write_network(output_path, relations, t) 

        self.outputs[0].description = path.basename(output_path)[:-4]
        
