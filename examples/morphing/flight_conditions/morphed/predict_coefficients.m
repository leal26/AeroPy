clc; clear all;
% Load predictor
addpath(genpath('\\coe-fs.engr.tamu.edu\Grads\leal26\Documents\GitHub\p3ga\kriging\dace'))

% Load Krigin for L/D ratio for morphing
% to predict L/D ratio predictor([alpha, V, H], NACA0012.performance_model)
load('\\coe-fs.engr.tamu.edu\Grads\leal26\Documents\GitHub\morphing_airfoils\Main\case_studies\comparisons\comparisons.mat')

% Code for lift coefficient Krigin
addpath(genpath('\\coe-fs.engr.tamu.edu\grads\leal26\Documents\GitHub\morphing_airfoils\Main\postprocessing'))
% [CL]=get_CL('NACA0012',alpha,V,H)

% Load angle of attack and velocity from file generated by opensky
load('./flight_conditions.mat')

% alpha = theta(:,:,1);
% V = theta(:,:,2);
% H = 10000*ones(size(alpha));

alpha = linspace(0,12,200);
v = linspace(20, 65, 200);
[alpha_grid, v_grid] = meshgrid(alpha, v);
H_grid = 10000*ones(size(alpha_grid));
alpha_grid = reshape(alpha_grid,[],1);
v_grid = reshape(v_grid,[],1);
H_grid = reshape(H_grid,[],1);

% NACA 0012
cl = get_CL('NACA0012',alpha_grid, v_grid,H_grid);
LD_ratio = predictor([alpha_grid, v_grid], NACA0012.performance_model);
fileID = fopen('morphed_NACA0012.txt','w');
formatSpec = '%f\t%f\t%f\t%f\n';
fprintf(fileID, formatSpec, [alpha_grid, v_grid, cl, LD_ratio].');
fclose(fileID);

% NACA 4415
cl = get_CL('NACA4415',alpha_grid, v_grid,H_grid);
LD_ratio = predictor([alpha_grid, v_grid], NACA4415.performance_model);
fileID = fopen('morphed_NACA4415.txt','w');
formatSpec = '%f\t%f\t%f\t%f\n';
fprintf(fileID, formatSpec, [alpha_grid, v_grid, cl, LD_ratio].');
fclose(fileID);

% NACA 641212
cl = get_CL('NACA641212',alpha_grid, v_grid,H_grid);
LD_ratio = predictor([alpha_grid, v_grid], NACA641212.performance_model);
fileID = fopen('morphed_NACA641212.txt','w');
formatSpec = '%f\t%f\t%f\t%f\n';
fprintf(fileID, formatSpec, [alpha_grid, v_grid, cl, LD_ratio].');
fclose(fileID);

% Glider
cl = get_CL('glider',alpha_grid, v_grid,H_grid);
LD_ratio = predictor([alpha_grid, v_grid], glider.performance_model);
fileID = fopen('morphed_glider.txt','w');
formatSpec = '%f\t%f\t%f\t%f\n';
fprintf(fileID, formatSpec, [alpha_grid, v_grid, cl, LD_ratio].');
fclose(fileID);
