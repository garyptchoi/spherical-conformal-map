function Nv = compute_vertex_normal(v,f)
% Compute the vertex normal of a mesh.
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nf = length(f);

f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);
e1 = v(f2,:) - v(f1,:);
e2 = v(f3,:) - v(f1,:);

% take cross product for each pair of successive edges
cross12 = cross(e1,e2);
area = abs(1/2*(cross12(:,1).^2+cross12(:,2).^2+cross12(:,3).^2).^(1/2));
% obtain face normal
Nf = [(1./(2*area)).*cross12(:,1), (1./(2*area)).*cross12(:,2), (1./(2*area)).*cross12(:,3)];


row = [f3; f1; f2];
col = [1:nf, 1:nf, 1:nf]';
val = [area; area; area]';
M = sparse(row,col,val,nv,nf);
% normalize
vertex_area_sum = sum(M,2);
[Mrow,Mcol,Mval] = find(M);
M = sparse(Mrow,Mcol,Mval./vertex_area_sum(Mrow),nv,nf);

% obtain vertex normal
Nv = M*Nf;