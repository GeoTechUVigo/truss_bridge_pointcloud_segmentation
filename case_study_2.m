clear all; clc; 

path_in='C:\Datos\Traballo\Codes\datos_publicar\case_study_2.las';
path_out = 'C:\Datos\Traballo\Codes\datos_publicar\case_study_2_segmented.las';
x_limits = [-16, 14.5];

width_faces_v_h_i = [0.6, 0.6, 0.6];
lateralBarWidth = 0.3;
verticalBarWidth = 0.2;
chordWidthZ = 0.5;

TrussSegmentation(path_in, path_out, width_faces_v_h_i, lateralBarWidth, verticalBarWidth, chordWidthZ, 'x_limits', x_limits, 'inner_horizontal', false, 'inside_board', true)
