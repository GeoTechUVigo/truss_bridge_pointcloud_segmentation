# Automated instance and semantic segmentation of point clouds of large metallic truss bridges with modelling purposes

Created by [Daniel Lamas Novoa](https://orcid.org/0000-0001-7275-183X), [Andrés Justo Dominguez](https://orcid.org/0000-0003-2072-4076), [Mario Soilán Rodríguez](https://orcid.org/0000-0001-6545-2225), [Manuel Cabaleiro Núñez](https://orcid.org/0000-0002-6948-1389) and [Belén Riveiro Rodríguez](https://orcid.org/0000-0002-1497-4370) from [GeoTech Group](https://geotech.webs.uvigo.es/en/), [CINTECX](http://cintecx.uvigo.es/gl/), [UVigo](https://www.uvigo.gal/).

## Overview
This repository contains the code of the work entitled [INSERTAR ENLACE ARTIGO].
Consists of a methodology for automatically segment point cloud truss bridges among 2 pillars and without them. The segmentation includes not only semantic but also instance information. The aim of the method is not segmenting all the points of the point cloud, but to ensure that each bar is segmented as an independent instance, avoiding that several bars are segmented as a single instance.
The method is a heuristic method which follows the following diagram.

![main_diagram_2](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/main_diagram_2.png)

In order to use the same functions, each face is reoriented as a vertical face. The steps applied to each face are shown in the following diagram.
![analysing_faces](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/analysing_faces.png)

Here is presented one of the case studies used to validate the methodology.
![bridge_segmentation](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/Images/bridge_segmentation.gift)


## Citation

## Licence
The code is released under WTFPL License (see [LICENSE](https://github.com/GeoTechUVigo/truss_bridge_pointcloud_segmentation/blob/main/LICENSE) file for details).

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
- path_out: path the LAS output file.
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
