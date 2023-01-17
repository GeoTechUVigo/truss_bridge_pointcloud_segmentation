%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
This file is part of the program Automated instance and semantic 
segmentation of point clouds of large metallic truss bridges with modelling
purposes.
The program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) 
any later version.
The program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details.
You should have received a copy of the GNU General Public License along 
with the program in COPYING. If not, see <https://www.gnu.org/licenses/>. 
%}

clear all; clc; 

path_in='C:\Datos\Traballo\Codes\datos_publicar\case_study_2.las';
path_out = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_2_segmented.las';

width_faces_v_h_i = [0.6, 0.6, 0.6];
lateralBarWidth = 0.3;
verticalBarWidth = 0.2;
chordWidthZ = 0.5;

TrussSegmentation(path_in, path_out, width_faces_v_h_i, lateralBarWidth, verticalBarWidth, chordWidthZ, 'inner_horizontal', false, 'inside_board', true)
