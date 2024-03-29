# [Instance and semantic segmentation of point clouds of large metallic truss bridges](https://doi.org/10.1016/j.autcon.2023.104865)

Created by [Daniel Lamas Novoa](https://orcid.org/0000-0001-7275-183X), [Andrés Justo Dominguez](https://orcid.org/0000-0003-2072-4076), [Mario Soilán Rodríguez](https://orcid.org/0000-0001-6545-2225), [Manuel Cabaleiro Núñez](https://orcid.org/0000-0002-6948-1389) and [Belén Riveiro Rodríguez](https://orcid.org/0000-0002-1497-4370) from [GeoTech Group](https://geotech.webs.uvigo.es/en/), [CINTECX](http://cintecx.uvigo.es/gl/), [UVigo](https://www.uvigo.gal/).

## Overview
This repository contains the code of the paper entitled [Instance and semantic segmentation of point clouds of large metallic truss bridges](https://doi.org/10.1016/j.autcon.2023.104865)

Consists of a methodology for automatically segment point cloud truss bridges among 2 pillars and without them. The segmentation includes not only semantic but also instance information. The aim of the method is not segmenting all the points of the point cloud, but to ensure that each bar is segmented as an independent instance, avoiding that several bars are segmented as a single instance.

The input required is a truss bridge point cloud LAS file. The output is the same file, with the semantic segmentation written in the ```Classification``` field and the instance segmentation in the ```PointSoourceId``` field.

The method is a heuristic method which follows the following diagram.

![main_diagram_2](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/main_diagram.png)

In order to use the same functions, each face is reoriented as a vertical face. The steps applied to each face are shown in the following diagram.
![analysing_faces](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/analysing_faces.png)

Here is presented one of the case studies used to validate the methodology.
![bridge_segmentation](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/bridge_segmentation.gif)


## Citation
If you find our work useful in your research, please consider citing:
```
@article{LAMAS2023104865,
title = {Instance and semantic segmentation of point clouds of large metallic truss bridges},
journal = {Automation in Construction},
volume = {151},
pages = {104865},
year = {2023},
issn = {0926-5805},
doi = {https://doi.org/10.1016/j.autcon.2023.104865},
url = {https://www.sciencedirect.com/science/article/pii/S0926580523001255},
author = {Daniel Lamas and Andrés Justo and Mario Soilán and Manuel Cabaleiro and Belén Riveiro},
keywords = {Point clouds, Truss bridge, Panoptic segmentation, Semantic segmentation, Instance segmentation, LiDAR, TLS},
abstract = {Several methods have been developed for the semantic segmentation of reinforced concrete bridges, however, there is a gap for truss bridges. Therefore, in this study a state-of-the-art methodology for the instance and semantic segmentation of point clouds of truss bridges for modelling purposes is presented, which, to the best of the authors' knowledge, is the first such methodology. This algorithm segments each truss element and classifies them as a chord, diagonal, vertical post, interior lateral brace, bottom lateral brace, or strut. The algorithm consists of a sequence of methods, including principal component analysis or clustering, that analyse each point and its neighbours in the point cloud. Case studies show that by adjusting only six manually measured parameters, the algorithm can automatically segment a truss bridge point cloud.}
}
```


## Licence
Automated instance and semantic segmentation of point clouds of large metallic truss bridges.

Copyright (C) 2023 GeoTECH Group <geotech@uvigo.gal>

Copyright (C) 2023 Daniel Lamas Novoa <daniel.lamas.novoa@uvigo.gal>

Copyright (C) 2023 Andrés Justo Domínguez <andres.justo.dominguez@uvigo.gal>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program in ![COPYING](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/COPYING.md). If not, see <https://www.gnu.org/licenses/>. 

## Data
The data is available in the [GeoTech Group repository](https://universidadevigo-my.sharepoint.com/:f:/g/personal/geotech_uvigo_gal/EoT3-ehKcexOs0yT2zS_LpABNX2Y-rswZvqBOB5cAgtt0Q), in the folder:

```
/DataSets/trus_bridge_pc/
```

## Requirements
The code runs in MATLAB2021b. The following Add-Ons are required:

- Circle fit 1.0.0.0
- Symbolic Math Toolbox 9.0
- Signal Processing Toolbox 8.7
- Statistics and Machine Learning Toolbox 12.2
- Image Processsing Toolbox 11.4
- Computer Vision Toolbox 10.1
 - Computer Vision Toolbox 10.1

## Usage
The main function is ![TrussSegmentation.m](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/TrussSegmentation.m)
In order to validate this methodology, 2 case studies are presented.
The data available for download used to test the methodology are already segmented. Evidently, the segmentation fields are not used during the segmnetation and they can be remove to carry out the test.

### Test case studies
Download the point cloud ```case_study_1_segmented.laz``` and the alignment ```case_study_1_alginment.csv``` for the 1st case study or the ```case_study_2_segmented.laz``` for the second. Notice that the 2nd does not have alignment file.

Unzip LAZ in LAS (for instance, using LAStools)

Open ![case_study_1.m](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/case_study_1.m) with MATLAB for the 1st or the ![case_study_2.m](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/case_study_2.m) for the 2nd case study.

Replace the value of the variables ```path_in``` and ```alignmentFile``` with your paths. Write your ```path_out``` value.

The parameters are already adjusted for this case study.

Run the code.

### Your on case study
Run the TrussSegmentation function adjusting the input parameters to your case study. The inputs of the function are:

- path_in: path the LAS input file.
- path_out: path the LAS output file. Semantic segmentation written in the ```Classification``` field and the instance segmentation in the ```PointSoourceId``` field.
- width_faces_v_h_i: 1x3 with the width of vertical, horizontal and inner faces.
- lateralBarWidth: diagonal and lateral brace width.
- verticalBarWidth: vertical bar width.
- chordWidthZ: chord width.
- grid: grid for segmentation. Default: 0.05.
- inside_board: True if the bridge board is not in the top of the structure.
- maxDistCentralFace: distance in the main direction of the bridge of 2 vertical members to be consider part of the same central face.
- numHorizontalFaces: number ot horizontal faces.
- numVerticalFaces: number of vertical faces.
- x_limits: not used.
- alignmentFile: path of the CSV alginment file. Not required.
- plot_vx: if True, the voxelised and segmented point cloud is shown.
- plot: if true, the point cloud segmented is shown.
- inner_vertical: True if threre are inner vertical bars.
- inner_horizontal: True if there are inner horizontal bars.
- path_out_csv: not used.
