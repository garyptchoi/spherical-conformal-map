function map = spherical_tutte_map(f, bigtri)
% Compute the spherical Tutte map using the approach in [1], with the cotangent Laplacian replaced by the Tutte Laplacian.
% Invoked only if the harmonic map fails due to very bad triangulations.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

if nargin < 2
    bigtri = 1;
end

nv = max(max(f)); 
nf = length(f);

% Construct the Tutte Laplacian
I = reshape(f',nf*3,1);
J = reshape(f(:,[2 3 1])',nf*3,1);
V = ones(nf*3,1)/2;
W = sparse([I;J],[J;I],[V;V]);
M = W + sparse(1:nv,1:nv,-diag(W)-(sum(W)'), nv, nv);

boundary = f(bigtri,1:3);
[mrow,mcol,mval] = find(M(boundary,:));
M = M - sparse(boundary(mrow), mcol, mval, nv, nv) + sparse(boundary, boundary, [1,1,1], nv, nv);
    
% set the boundary condition for big triangle
b = zeros(nv,1);
b(boundary) = exp(1i*(2*pi*(0:2)/length(boundary)));

% solve the Laplace equation to obtain a Tutte map
z = M \ b;
z = z-mean(z);

% inverse stereographic projection
S = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), (-1+abs(z).^2)./(1+abs(z).^2)];

%% Find optimal big triangle size
w = complex(S(:,1)./(1+S(:,3)), S(:,2)./(1+S(:,3)));

% find the index of the southernmost triangle
[~, index] = sort(abs(z(f(:,1)))+abs(z(f(:,2)))+abs(z(f(:,3))));
inner = index(1);
if inner == bigtri
    inner = index(2);
end

% Compute the size of the northern most and the southern most triangles 
NorthTriSide = (abs(z(f(bigtri,1))-z(f(bigtri,2))) + ...
    abs(z(f(bigtri,2))-z(f(bigtri,3))) + ...
    abs(z(f(bigtri,3))-z(f(bigtri,1))))/3;

SouthTriSide = (abs(w(f(inner,1))-w(f(inner,2))) + ...
    abs(w(f(inner,2))-w(f(inner,3))) + ...
    abs(w(f(inner,3))-w(f(inner,1))))/3;

% rescale to get the best distribution
z = z*(sqrt(NorthTriSide*SouthTriSide))/(NorthTriSide); 

% inverse stereographic projection
map = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), (-1+abs(z).^2)./(1+abs(z).^2)];

