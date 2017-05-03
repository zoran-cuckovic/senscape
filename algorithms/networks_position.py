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
                                        ParameterBoolean,
                                        ParameterString,
                                        ParameterSelection )
from processing.core.outputs import  OutputRaster ###OutputVector
from processing.tools import dataobjects


class NetworkPositions(GeoAlgorithm):
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


    INPUT_DEM = 'INPUT_DEM'
    INPUT_TOPONET_RIDGE = 'INPUT_TOPONET_RIDGE'
    INPUT_TOPONET_VALLEY = 'INPUT_TOPONET_VALLEY'

    OUTPUT_RASTER = 'OUTPUT_RASTER'

    
    def defineCharacteristics(self):
        """Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # The name that the user will see in the toolbox
        self.name = 'Topographic position'

        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Topographic networks'

        # We add the input vector layer. It can have any kind of geometry
        # It is a mandatory (not optional) one, hence the False argument
	
		
        self.addParameter(ParameterRaster(self.INPUT_DEM,
            self.tr('Elevation model (DEM)'),  False))
        
        self.addParameter(ParameterRaster(self.INPUT_TOPONET_RIDGE,
            self.tr('Ridge network'),  False))
        
        self.addParameter(ParameterRaster(self.INPUT_TOPONET_VALLEY,
            self.tr('Valley network'),  False))
        
        self.addOutput(OutputRaster(self.OUTPUT_RASTER,
                                    self.tr('Output layer')))


        
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
        
	
        import gdal
        import numpy as np
        from .modules import doNetworks as nt


        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        Raster_path = self.getParameterValue(self.INPUT_DEM)
        Ridges_path= self.getParameterValue(self.INPUT_TOPONET_RIDGE)
        Valleys_path= self.getParameterValue(self.INPUT_TOPONET_VALLEY)
        
        Output_raster = self.getOutputValue(self.OUTPUT_RASTER)

        d= gdal.Open(Raster_path)
        r=d.ReadAsArray()

        
        for s in [1, 2]:

            if s ==1: d2 = gdal.Open(Ridges_path)
            else: gdal.Open(Valleys_path)


            ids, indices, steps = d2.ReadAsArray()

            ids = None 

            y_parent, x_parent = np.unravel_index(indices, indices.shape)

            slopes = nt.slope( r, x_parent, y_parent)
            
            degrees = nt.degrees(x_parent, y_parent)

  
            if s == 1 :
                out = nt.accum_slope(slopes, x_parent, y_parent, steps, degrees)
            else:
                # attention! Valley slopes are negative, so we go += to subtract
                out += nt.accum_slope(slopes, x_parent, y_parent, steps, degrees)

        out= nt.sinks(x_parent, y_parent, ids, borders=True)

        
        
        #out = nt.accum_min(r, x_parent, y_parent, steps)
# DOES NOT WORK PROPERLY ??? does nto catch child_ids for output

        driver = gdal.GetDriverByName('GTiff')
##
##   
        c = driver.CreateCopy(Output_raster, d,0) #WHAT IS 0 ? : "strict or not" default =1
        c.GetRasterBand(1).WriteArray(out)
##

        
        
                        
        
