% Demo of extending our proposed spherical mapping algorithms for
% minimizing other distortions, such as local area distortion
%
% In case you want to focus on spherical conformal map, read demo.m instead
%
% If you use this code in your own work, please cite the following papers:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% (For mobius_area_correction_spherical)
% [2] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
%     "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
%     SIAM Journal on Imaging Sciences, 2020.

% Copyright (c) 2013-2020, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

addpath('mfile')
addpath('extension') % contain the codes for area-preserving map

%% Example 1: David
load('david.mat')
plot_mesh(v,f); view([-130 0])

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f); title('Spherical conformal map')

%% Extension 1: a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f); title('Spherical area-preserving map')

%% Extension 2: an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f); title('Iterative spherical area-preserving map')

%% Extension 3: our linear method for spherical conformal map together with a Mobius area correction (Choi et al., SIAM J. Imaging Sci. 2020)
map = spherical_conformal_map(v,f);
map = mobius_area_correction_spherical(v,f,map);
plot_mesh(map,f); title('Spherical conformal map with Mobius area correction')

%% evaluate the angle distortion (please run any of the above methods and then run this part for the evaluation)
d = angle_distortion(v,f,map);
a = area_distortion(v,f,map);

fprintf('Mean(angle distortion) = %.4f\n',mean(abs(d)));
fprintf('SD(angle distortion) = %.4f\n',std(abs(d)));
fprintf('Mean(area distortion) = %.4f\n',mean(abs(a)));
fprintf('SD(area distortion) = %.4f\n',std(abs(a)));




%% Example 2: Lion
load('lion.mat')
plot_mesh(v,f,mean_curv);

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Spherical conformal map')

%% Extension 1: a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Spherical area-preserving map')

%% Extension 2: an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Iterative spherical area-preserving map')

%% Extension 3: our linear method for spherical conformal map together with a Mobius area correction
map = spherical_conformal_map(v,f);
map = mobius_area_correction_spherical(v,f,map);
plot_mesh(map,f,mean_curv); view([-70 0]); title('Spherical conformal map with Mobius area correction')




%% Example 3: Brain
load('brain.mat')
plot_mesh(v,f,mean_curv); view([90 0]); 

%% our linear method for spherical conformal map
map = spherical_conformal_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Spherical conformal map')

%% Extension 1: a simple variation of our linear method for a more area-preserving map
map = spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Spherical area-preserving map')

%% Extension 2: an iterative scheme for further reducing the area distortion
map = iterative_spherical_area_preserving_map(v,f);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Iterative spherical area-preserving map')

%% Extension 3: our linear method for spherical conformal map together with a Mobius area correction
map = spherical_conformal_map(v,f);
map = mobius_area_correction_spherical(v,f,map);
plot_mesh(map,f,mean_curv); view([-30 0]); title('Spherical conformal map with Mobius area correction')