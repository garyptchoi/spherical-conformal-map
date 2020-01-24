% Demo of extending our proposed spherical mapping algorithms for
% minimizing other distortions, such as local area distortion
%
% In case you want to focus on spherical conformal map, read demo.m instead
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

addpath('mfile')
addpath('extension') % contain the codes for area-preserving map

%% Example 1: David
load('david.mat')
plot_mesh(v,f); view([-130 0])

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f); title('Spherical conformal map')

%% a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f); title('Spherical area-preserving map')

%% an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f); title('Iterative spherical area-preserving map')

%%
% evaluate the angle distortion
% d = angle_distortion(v,f,map);
a = area_distortion(v,f,map);
[ mean(abs(a)), std(abs(a))]
% [mean(abs(d)), std(abs(d)), mean(abs(a)), std(abs(a))]


%% Example 2: Lion
load('lion.mat')
plot_mesh(v,f,mean_curv);

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Spherical conformal map')

%% a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Spherical area-preserving map')

%% an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Iterative spherical area-preserving map')



%% Example 3: Brain
load('brain.mat')
plot_mesh(v,f,mean_curv); view([90 0]); 

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Spherical conformal map')

%% a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Spherical area-preserving map')

%% an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Iterative spherical area-preserving map')
