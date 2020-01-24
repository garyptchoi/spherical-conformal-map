function map = linear_beltrami_solver(v,f,mu,landmark,target)
% Linear Beltrami solver.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

af = (1-2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);
bf = -2*imag(mu)./(1.0-abs(mu).^2);
gf = (1+2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);

f0 = f(:,1); f1 = f(:,2); f2 = f(:,3);

uxv0 = v(f1,2) - v(f2,2);
uyv0 = v(f2,1) - v(f1,1);
uxv1 = v(f2,2) - v(f0,2);
uyv1 = v(f0,1) - v(f2,1); 
uxv2 = v(f0,2) - v(f1,2);
uyv2 = v(f1,1) - v(f0,1);

l = [sqrt(sum(uxv0.^2 + uyv0.^2,2)), sqrt(sum(uxv1.^2 + uyv1.^2,2)), sqrt(sum(uxv2.^2 + uyv2.^2,2))];
s = sum(l,2)*0.5;

area = sqrt(s.*(s-l(:,1)).*(s-l(:,2)).*(s-l(:,3)));

v00 = (af.*uxv0.*uxv0 + 2*bf.*uxv0.*uyv0 + gf.*uyv0.*uyv0)./area;
v11 = (af.*uxv1.*uxv1 + 2*bf.*uxv1.*uyv1 + gf.*uyv1.*uyv1)./area;
v22 = (af.*uxv2.*uxv2 + 2*bf.*uxv2.*uyv2 + gf.*uyv2.*uyv2)./area;
v01 = (af.*uxv1.*uxv0 + bf.*uxv1.*uyv0 + bf.*uxv0.*uyv1 + gf.*uyv1.*uyv0)./area;
v12 = (af.*uxv2.*uxv1 + bf.*uxv2.*uyv1 + bf.*uxv1.*uyv2 + gf.*uyv2.*uyv1)./area;
v20 = (af.*uxv0.*uxv2 + bf.*uxv0.*uyv2 + bf.*uxv2.*uyv0 + gf.*uyv0.*uyv2)./area;

I = [f0;f1;f2;f0;f1;f1;f2;f2;f0];
J = [f0;f1;f2;f1;f0;f2;f1;f0;f2];
V = [v00;v11;v22;v01;v01;v12;v12;v20;v20]/2;
A = sparse(I,J,-V);

targetc = target(:,1) + 1i*target(:,2);
b = -A(:,landmark)*targetc;
b(landmark) = targetc;
A(landmark,:) = 0; A(:,landmark) = 0;
A = A + sparse(landmark,landmark,ones(length(landmark),1), size(A,1), size(A,2));
map = A\b;
map = [real(map),imag(map)];

end
