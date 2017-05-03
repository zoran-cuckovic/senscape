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
from processing.core.outputs import  OutputRaster
from processing.tools import dataobjects


class NetworkMinAccum(GeoAlgorithm):
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


    INPUT_DEM = 'INPUT_RASTER'
    INPUT_TOPONET= 'INPUT_TOPONET'
    DIRECTION = 'DIRECTION'

    OUTPUT_RASTER = 'OUTPUT_RASTER'

    
    def defineCharacteristics(self):
        """Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # The name that the user will see in the toolbox
        self.name = 'Accum minimum vs max'

        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Topographic networks'

        # We add the input vector layer. It can have any kind of geometry
        # It is a mandatory (not optional) one, hence the False argument
	
		
        self.addParameter(ParameterRaster(self.INPUT_DEM,
            self.tr('Elevation model (DEM)'),  False))

        self.addParameter(ParameterRaster(self.INPUT_TOPONET,
            self.tr('Topographic network'),  False))

        self.addParameter(ParameterSelection(self.DIRECTION,
                                            self.tr('Direction'),
                ['Lowest', 'Highest'], 0))

        self.addOutput(OutputRaster(self.OUTPUT_RASTER,
                                    self.tr('Output raster')))


        
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
        from .modules import TopoGraph as tg


        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        Raster_path = self.getParameterValue(self.INPUT_DEM)
        Network_path= self.getParameterValue(self.INPUT_TOPONET)
        
        Direction = self.getParameterValue(self.DIRECTION)
        
        Output_raster = self.getOutputValue(self.OUTPUT_RASTER)

        d= gdal.Open(Raster_path)
        r=d.ReadAsArray()


        if Direction: r*=-1


        n= gdal.Open(Network_path)
        ids, indices, steps=n.ReadAsArray()

        y_parent, x_parent = np.unravel_index(indices, indices.shape)
      
        print "accumulation"
        mins = nt.accum_min (r,  x_parent, y_parent, steps)

        if Direction: mins *= -1 # to restore normal values...
        
        driver = gdal.GetDriverByName('GTiff')

        # FORMAT : should be flexible: int for ranks etc
        outDs = driver.Create(Output_raster,
                              d.RasterXSize, d.RasterYSize,1,
                              gdal.GDT_Float32)       
        outDs.SetGeoTransform(d.GetGeoTransform())
        outDs.SetProjection(d.GetProjection())
       
        outDs.GetRasterBand(1).WriteArray(mins)
    
        
     
        
