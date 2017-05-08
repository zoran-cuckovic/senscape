# Senspace

Senspace is a toolbox for modelling human perception of topography for QGIS, based on [QGIS viewshed module](https://github.com/zoran-cuckovic/QGIS-visibility-analysis). It is a subclass of QGIS processing framework and can be readily used in custom models or scripts.

The plugin is in experimental phase.

**To install:** see manual instalation of [QGIS viewshed module](https://github.com/zoran-cuckovic/QGIS-visibility-analysis). It should be placed in a folder named "Senspace", and enabled under QGIS processing options (besides being enabled as a plugin in QGIS plugin manager).

The plugin is being developed with help by Alexander Bruy. 

## Instructions

The analysis is made in two steps:
- 1 Vector layer containing observer points is first cleaned and prepared using "Create viewpoints" module
- 2 Viewshed analysis is made with the "Viewshed" module, using the created points
- 2b For intervisibility analysis additional target points need to be created using the same module as in (1). Target height is specified in the corresponding text-box or by choosing a table field. 