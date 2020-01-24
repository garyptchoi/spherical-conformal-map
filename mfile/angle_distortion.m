function distortion = angle_distortion(v,f,map)

% Calculate and visualize the angle difference (angle_map - angle_v)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% map: nv x 3 vertex coordinates of the mapping result
%
% Output:
% distortion: 3*nf x 1 angle differences
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nv2 = length(map);

if nv ~= nv2
    error('Error: The two meshes are of different size.');
end

if size(v,2) == 1
    v = [real(v),imag(v),zeros(length(v),1)];
elseif size(v,2) == 2
    v = [v,zeros(length(v),1)];
end

if size(map,2) == 1
    map = [real(map),imag(map),zeros(length(map),1)];
elseif size(map,2) == 2
    map = [map,zeros(length(map),1)];
end

f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

% calculate angles on v
a3=[v(f1,1)-v(f3,1), v(f1,2)-v(f3,2), v(f1,3)-v(f3,3)];
b3=[v(f2,1)-v(f3,1), v(f2,2)-v(f3,2), v(f2,3)-v(f3,3)];
a1=[v(f2,1)-v(f1,1), v(f2,2)-v(f1,2), v(f2,3)-v(f1,3)];
b1=[v(f3,1)-v(f1,1), v(f3,2)-v(f1,2), v(f3,3)-v(f1,3)];
a2=[v(f3,1)-v(f2,1), v(f3,2)-v(f2,2), v(f3,3)-v(f2,3)];
b2=[v(f1,1)-v(f2,1), v(f1,2)-v(f2,2), v(f1,3)-v(f2,3)];
vcos1=(a1(:,1).*b1(:,1)+a1(:,2).*b1(:,2)+a1(:,3).*b1(:,3))./ ...
    (sqrt(a1(:,1).^2+a1(:,2).^2+a1(:,3).^2).*sqrt(b1(:,1).^2+b1(:,2).^2+b1(:,3).^2));
vcos2=(a2(:,1).*b2(:,1)+a2(:,2).*b2(:,2)+a2(:,3).*b2(:,3))./...
    (sqrt(a2(:,1).^2+a2(:,2).^2+a2(:,3).^2).*sqrt(b2(:,1).^2+b2(:,2).^2+b2(:,3).^2));
vcos3=(a3(:,1).*b3(:,1)+a3(:,2).*b3(:,2)+a3(:,3).*b3(:,3))./...
    (sqrt(a3(:,1).^2+a3(:,2).^2+a3(:,3).^2).*sqrt(b3(:,1).^2+b3(:,2).^2+b3(:,3).^2));
    
% calculate angles on map
c3=[map(f1,1)-map(f3,1), map(f1,2)-map(f3,2), map(f1,3)-map(f3,3)];
d3=[map(f2,1)-map(f3,1), map(f2,2)-map(f3,2), map(f2,3)-map(f3,3)];
c1=[map(f2,1)-map(f1,1), map(f2,2)-map(f1,2), map(f2,3)-map(f1,3)];
d1=[map(f3,1)-map(f1,1), map(f3,2)-map(f1,2), map(f3,3)-map(f1,3)];
c2=[map(f3,1)-map(f2,1), map(f3,2)-map(f2,2), map(f3,3)-map(f2,3)];
d2=[map(f1,1)-map(f2,1), map(f1,2)-map(f2,2), map(f1,3)-map(f2,3)];
mapcos1=(c1(:,1).*d1(:,1)+c1(:,2).*d1(:,2)+c1(:,3).*d1(:,3))./...
    (sqrt(c1(:,1).^2+c1(:,2).^2+c1(:,3).^2).*sqrt(d1(:,1).^2+d1(:,2).^2+d1(:,3).^2));
mapcos2=(c2(:,1).*d2(:,1)+c2(:,2).*d2(:,2)+c2(:,3).*d2(:,3))./...
    (sqrt(c2(:,1).^2+c2(:,2).^2+c2(:,3).^2).*sqrt(d2(:,1).^2+d2(:,2).^2+d2(:,3).^2));
mapcos3=(c3(:,1).*d3(:,1)+c3(:,2).*d3(:,2)+c3(:,3).*d3(:,3))./...
    (sqrt(c3(:,1).^2+c3(:,2).^2+c3(:,3).^2).*sqrt(d3(:,1).^2+d3(:,2).^2+d3(:,3).^2));
    
% calculate the angle difference
distortion = (acos([mapcos1;mapcos2;mapcos3])-acos([vcos1;vcos2;vcos3]))*180/pi;

% histogram
figure;
hist(distortion,-180:1:180);
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);
