function E = evaluate_energy(v, M)
% Evaluate the spring energy with the weighted Laplacian matrix M.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

[row, col, value] = find(M);
diagonal = find(row>=col);
row(diagonal) = [];
col(diagonal) = [];
value(diagonal) = [];
d = v(row,1:3)-v(col,1:3);
E = sum(value.*(d(:,1).^2 + d(:,2).^2 + d(:,3).^2));
E = abs(sum(E))/2;