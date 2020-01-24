% A linear method for computing spherical conformal map of genus-0 closed triangle meshes
%
% Main program:
% map = spherical_conformal_map(v,f)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
%
% Output:
% map: nv x 3 vertex coordinates of the spherical conformal map
% 
% Remark:
% See demo_extension.m to understand how the method can be extended for
% achieving other effects, e.g. spherical area-preserving map
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

addpath('mfile')

%% Example 1: David
load('david.mat')

plot_mesh(v,f); view([-130 0])

map = spherical_conformal_map(v,f);

plot_mesh(map,f); 

% evaluate the angle distortion
angle_distortion(v,f,map);

%% Example 2: Chinese lion
load('lion.mat')
% plot_mesh(v,f);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(v,f,mean_curv);

map = spherical_conformal_map(v,f);

% plot_mesh(map,f); view([-70 0])
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); view([-70 0])

% evaluate the angle distortion
angle_distortion(v,f,map);

%% Example 3: Brain
load('brain.mat')
% plot_mesh(v,f); view([90 0])
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(v,f,mean_curv); view([90 0]); 

map = spherical_conformal_map(v,f);

% plot_mesh(map,f); view([-30 0]);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); view([-30 0]);

% evaluate the angle distortion
angle_distortion(v,f,map);