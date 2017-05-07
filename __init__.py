# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Senspace
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
 This script initializes the plugin, making it known to QGIS.
"""

__author__ = 'Zoran Čučković'
__date__ = '2017-05-01'
__copyright__ = '(C) 2017 by Zoran Čučković'

from senspace_plugin import SenspacePlugin

# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load Senspace class from file Senspace.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    
    return SenspacePlugin()
