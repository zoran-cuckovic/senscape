# -*- coding: utf-8 -*-

"""
/***************************************************************************
 Senscape
                                 A QGIS plugin
 Description to fill
                              -------------------
        begin                : 2017-05-01
        copyright            : (C) 2017 by Zoran Čučković
        email                : N/A
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

__author__ = 'Zoran Čučković'
__date__ = '2017-05-01'
__copyright__ = '(C) 2017 by Zoran Čučković'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from processing.core.AlgorithmProvider import AlgorithmProvider
from processing.core.ProcessingConfig import Setting, ProcessingConfig


from .algorithms.viewshed_raster  import  Viewshed
from .algorithms.viewshed_points import  ViewshedPoints
from .algorithms.viewshed_intervisibility import  Intervisibility

#from .algorithms.networks_position import  NetworkPositions
#from .algorithms.networks_rank import NetworkRanks
#from .algorithms.networks_create import NetworkCreate
#from .algorithms.networks_min_accum import NetworkMinAccum
#from .algorithms.networks_vector import NetworkVector



from PyQt4.QtGui import QIcon
from os import path


class SenscapeProvider(AlgorithmProvider):

    MY_DUMMY_SETTING = 'MY_DUMMY_SETTING'

    def __init__(self):
        AlgorithmProvider.__init__(self)

        # Deactivate provider by default
        self.activate = True

        # Load algorithms
        self.alglist = [ Intervisibility(), Viewshed(), ViewshedPoints()  ]
        """
        TODO
			NetworkPositions(), 
			NetworkCreate(), NetworkRanks(),
			NetworkMinAccum(), NetworkVector()]
        """
        for alg in self.alglist:
            alg.provider = self

    def initializeSettings(self):
        """In this method we add settings needed to configure our
        provider.

        Do not forget to call the parent method, since it takes care
        or automatically adding a setting for activating or
        deactivating the algorithms in the provider.
        """
        AlgorithmProvider.initializeSettings(self)
        ProcessingConfig.addSetting(Setting('Example algorithms',
            SenscapeProvider.MY_DUMMY_SETTING,
            'Example setting', 'Default value'))

    def unload(self):
        """Setting should be removed here, so they do not appear anymore
        when the plugin is unloaded.
        """
        AlgorithmProvider.unload(self)
        ProcessingConfig.removeSetting(
            SenscapeProvider.MY_DUMMY_SETTING)

    def getName(self):
	#def id(self):
	
        """This is the name that will NOT!! appear on the toolbox group.

        It is also used to create the command line name of all the
        algorithms from this provider.
        """
        return 'senscape'

    def getDescription(self):
	#def name (self):
        """This is the provired full name.
        """
        return 'Senscape'

    def getIcon(self):

        """We return the default icon.
        
	"""
        
        return QIcon(path.dirname(__file__) + '/icon.png')

        #return AlgorithmProvider.getIcon(self)


    def _loadAlgorithms(self):
        """Here we fill the list of algorithms in self.algs.

        This method is called whenever the list of algorithms should
        be updated. If the list of algorithms can change (for instance,
        if it contains algorithms from user-defined scripts and a new
        script might have been added), you should create the list again
        here.

        In this case, since the list is always the same, we assign from
        the pre-made list. This assignment has to be done in this method
        even if the list does not change, since the self.algs list is
        cleared before calling this method.
        """
        self.algs = self.alglist
