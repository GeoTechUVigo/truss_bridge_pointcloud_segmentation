# Automated instance and semantic segmentation of point clouds of large metallic truss bridges with modelling purposes

Created by [Daniel Lamas Novoa](https://orcid.org/0000-0001-7275-183X), [Andrés Justo Dominguez](https://orcid.org/0000-0003-2072-4076), [Mario Soilán Rodríguez](https://orcid.org/0000-0001-6545-2225), [Manuel Cabaleiro Núñez](https://orcid.org/0000-0002-6948-1389) and [Belén Riveiro Rodríguez](https://orcid.org/0000-0002-1497-4370) from [GeoTech Group](https://geotech.webs.uvigo.es/en/), [CINTECX](http://cintecx.uvigo.es/gl/), [UVigo](https://www.uvigo.gal/).

## Overview
This repository contains the code of the work entitled [INSERTAR ENLACE ARTIGO].
Consists of a methodology for automatically segment point cloud truss bridges. The segmentation includes not only semantic but also instance information. The aim of the method is not segmenting all the points of the point cloud, but to ensure that each bar is segmented as an independent instance, avoiding that several bars are segmented as a single instance.
The method is a heuristic method which follows the following diagram.

![Model](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/tree/main/Images/main_diagram_2.png)

In order to use the same functions, each face is reoriented as a vertical face. The steps applied to each face are shown in the following diagram.
![alt_text](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/tree/main/Images/analysing_faces.png)

Here is presented one of the case studies used to validate the methodology.
![alt_text](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/tree/main/Images/bridge_pc.png)
![alt_text](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/tree/main/Images/bridge_segmented.png)


## Citation


## Licence


## Usage


