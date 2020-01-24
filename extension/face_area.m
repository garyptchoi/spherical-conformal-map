function fa = face_area(f,v)
% Compute the area of every face of a triangle mesh.
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

v12 = v(f(:,2),:) - v(f(:,1),:);
v23 = v(f(:,3),:) - v(f(:,2),:);
v31 = v(f(:,1),:) - v(f(:,3),:);

a = sqrt(dot(v12,v12,2));
b = sqrt(dot(v23,v23,2));
c = sqrt(dot(v31,v31,2));

s = (a+b+c)/2;
fa = sqrt(s.*(s-a).*(s-b).*(s-c)); 
