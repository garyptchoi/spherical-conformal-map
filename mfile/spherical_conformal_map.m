function map = spherical_conformal_map(v,f)
% A linear method for computing spherical conformal map of a genus-0 closed surface
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
%
% Output:
% map: nv x 3 vertex coordinates of the spherical conformal parameterization
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%% Check whether the input mesh is genus-0
if length(v)-3*length(f)/2+length(f) ~= 2
    error('The mesh is not a genus-0 closed surface.\n');
end

%% Find the most regular triangle as the "big triangle"
temp = v(reshape(f',1,length(f)*3),1:3);
e1 = sqrt(sum((temp(2:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e2 = sqrt(sum((temp(1:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e3 = sqrt(sum((temp(1:3:end,1:3) - temp(2:3:end,1:3))'.^2))';
regularity = abs(e1./(e1+e2+e3)-1/3)+...
    abs(e2./(e1+e2+e3)-1/3)+abs(e3./(e1+e2+e3)-1/3);
[~,bigtri] = min(regularity);

% In case the spherical parameterization result is not evenly distributed,
% try to change bigtri to the id of some other triangles with good quality

%% North pole step: Compute spherical map by solving laplace equation on a big triangle
nv = size(v,1); 
M = cotangent_laplacian(v,f);

p1 = f(bigtri,1);
p2 = f(bigtri,2);
p3 = f(bigtri,3);

fixed = [p1,p2,p3];
[mrow,mcol,mval] = find(M(fixed,:));
M = M - sparse(fixed(mrow),mcol,mval,nv,nv) + sparse(fixed,fixed,[1,1,1],nv,nv);

% set the boundary condition for big triangle
x1 = 0; y1 = 0; x2 = 1; y2 = 0; % arbitrarily set the two points
a = v(p2,1:3) - v(p1,1:3);
b = v(p3,1:3) - v(p1,1:3);
sin1 = (norm(cross(a,b),2))/(norm(a,2)*norm(b,2));
ori_h = norm(b,2)*sin1;
ratio = norm([x1-x2,y1-y2],2)/norm(a,2);
y3 = ori_h*ratio; % compute the coordinates of the third vertex
x3 = sqrt(norm(b,2)^2*ratio^2-y3^2);

% Solve the Laplace equation to obtain a harmonic map
c = zeros(nv,1); c(p1) = x1; c(p2) = x2; c(p3) = x3;
d = zeros(nv,1); d(p1) = y1; d(p2) = y2; d(p3) = y3;
z = M \ complex(c,d);
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
S = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), (-1+abs(z).^2)./(1+abs(z).^2)];

if sum(sum(isnan(S))) ~= 0
    % if harmonic map fails due to very bad triangulations, use tutte map
    S = spherical_tutte_map(f,bigtri);
end

%% South pole step
[~,I] = sort(S(:,3));

% number of points near the south pole to be fixed  
% simply set it to be 1/10 of the total number of vertices (can be changed)
% In case the spherical parameterization is not good, change 10 to
% something smaller (e.g. 2)
fixnum = max(round(length(v)/10),3);
fixed = I(1:min(length(v),fixnum)); 

% south pole stereographic projection
P = [S(:,1)./(1+S(:,3)), S(:,2)./(1+S(:,3))]; 

% compute the Beltrami coefficient
mu = beltrami_coefficient(P, f, v); 

% compose the map with another quasi-conformal map to cancel the distortion
map = linear_beltrami_solver(P,f,mu,fixed,P(fixed,:)); 

if sum(sum(isnan(map))) ~= 0
    % if the result has NaN entries, then most probably the number of
    % boundary constraints is not large enough
    
    % increase the number of boundary constrains and run again
    fixnum = fixnum*5; % again, this number can be changed
    fixed = I(1:min(length(v),fixnum)); 
    map = linear_beltrami_solver(P,f,mu,fixed,P(fixed,:)); 
    
    if sum(sum(isnan(map))) ~= 0
        map = P; % use the old result
    end
end

z = complex(map(:,1),map(:,2));

% inverse south pole stereographic projection
map = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), -(abs(z).^2-1)./(1+abs(z).^2)];

end