%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
This file is part of the program Automated instance and semantic 
segmentation of point clouds of large metallic truss bridges with modelling
purposes.
The program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version.
The program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details.
You should have received a copy of the GNU General Public License along 
with the program in COPYING. If not, see <https://www.gnu.org/licenses/>. 
%}

clear all; clc; 

path_in='C:\Datos\Traballo\Codes\datos_publicar\case_study_1.las';
path_out = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_1_segmented.las';
alignmentFile = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_1_alignment.csv';

width_faces_v_h_i = [0.9, 0.6, 0.9];
lateralBarWidth = 0.1;
verticalBarWidth = 0.3;
chordWidthZ = 1;

TrussSegmentation(path_in, path_out, width_faces_v_h_i, lateralBarWidth, verticalBarWidth, chordWidthZ, 'alignmentFile', alignmentFile, 'inner_vertical', false)