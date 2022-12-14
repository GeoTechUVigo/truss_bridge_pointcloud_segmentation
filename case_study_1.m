clear all; clc; 

path_in='C:\Datos\Traballo\Codes\datos_publicar\case_study_1.las';
path_out = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_1_segmented.las';
alignmentFile = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_1_alignment.csv';
x_limits = [-43, 21.5];

width_faces_v_h_i = [0.9, 0.6, 0.9];
lateralBarWidth = 0.1;
verticalBarWidth = 0.3;
chordWidthZ = 1;

TrussSegmentation(path_in, path_out, width_faces_v_h_i, lateralBarWidth, verticalBarWidth, chordWidthZ, 'x_limits', x_limits, 'alignmentFile', alignmentFile, 'inner_vertical', false)