function [map_mobius,x] =  mobius_area_correction_spherical(v,f,map)

% Find an optimal Mobius transformation for reducing the area distortion of a spherical conformal parameterization using the method in [1].
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
% map: nv x 3 vertex coordinates of the spherical conformal parameterization
% 
% Output:
% map_mobius: nv x 3 vertex coordinates of the updated spherical conformal parameterization
% x: the optimal parameters for the Mobius transformation, where
%    f(z) = \frac{az+b}{cz+d}
%         = ((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i))
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
%     "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
%     SIAM Journal on Imaging Sciences, 2020.
%
% Copyright (c) 2019-2020, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%%
% Compute the area with normalization
area_v = face_area(f,v); area_v = area_v/sum(area_v);

% Project the sphere onto the plane
p = stereographic(map);
z = complex(p(:,1),p(:,2));


% Function for calculating the area after the Mobius transformation 
area_map = @(x) face_area(f,stereographic([real(((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i))),...
    imag(((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i)))]))/...
    sum(face_area(f,stereographic([real(((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i))),...
    imag(((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i)))])));

% objective function: mean(abs(log(area_map./area_v)))
d_area = @(x) finitemean(abs(log(area_map(x)./area_v)));

% Optimization setup
x0 = [1,0,0,0,0,0,1,0]; % initial guess
lb = [-1,-1,-1,-1,-1,-1,-1,-1]*100; % lower bound for the parameters
ub = [1,1,1,1,1,1,1,1]*100; % upper bound for the parameters
options = optimoptions('fmincon','Display','iter');

% Optimization (may further supply gradients for better result, not yet implemented)
x = fmincon(d_area,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = ((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i));
map_mobius = stereographic([real(fz), imag(fz)]);

end

function fa = face_area(f,v)
% Compute the area of every face of a triangle mesh.
v12 = v(f(:,2),:) - v(f(:,1),:);
v23 = v(f(:,3),:) - v(f(:,2),:);
v31 = v(f(:,1),:) - v(f(:,3),:);

a = sqrt(dot(v12,v12,2));
b = sqrt(dot(v23,v23,2));
c = sqrt(dot(v31,v31,2));

s = (a+b+c)/2;
fa = sqrt(s.*(s-a).*(s-b).*(s-c)); 
end

function m = finitemean(A)
% for avoiding the Inf values caused by division by a very small area
    m = mean(A(isfinite(A)));
end

function v = stereographic(u)
% STEREOGRAPHIC  Stereographic projection.
%   v = STEREOGRAPHIC(u), for N-by-2 matrix, projects points in plane to sphere
%                       ; for N-by-3 matrix, projects points on sphere to plane
    if size(u, 2) == 1
      u = [real(u), imag(u)];
    end
    x = u(:,1);
    y = u(:,2);
    if size(u,2) < 3
      z = 1 + x.^2 + y.^2;
      v = [2*x ./ z, 2*y ./ z, (-1 + x.^2 + y.^2) ./ z];
    else
      z = u(:,3);
      v = [x ./ (1-z), y ./ (1-z)];
    end
end