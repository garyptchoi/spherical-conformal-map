function M = mesh_laplacian(v,f,weight,normalization)
% Compute the mesh Laplacian with different weights.
% weight: authalic, Tutte, squared_inverse_distance, inverse_distance
% normalization: normalize the weights (default: 0)
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

if ~exist('normalization','var')
    normalization = 0; % no normalization
end

nv = length(v);

if strcmp(weight,'authalic')
    % The locally authalic weight
    % This weight is related to the Chi energy (Desbrun et al., 
    % Eurographics 2002), which measures local area distortion.

    f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

    l1 = sqrt(sum((v(f2,:) - v(f3,:)).^2,2));
    l2 = sqrt(sum((v(f3,:) - v(f1,:)).^2,2));
    l3 = sqrt(sum((v(f1,:) - v(f2,:)).^2,2));

    s = (l1 + l2 + l3)*0.5;
    area = sqrt(s.*(s-l1).*(s-l2).*(s-l3));

    w12 = (l1.^2 + l2.^2 - l3.^2)./area/4;
    w23 = (l2.^2 + l3.^2 - l1.^2)./area/4; 
    w31 = (l1.^2 + l3.^2 - l2.^2)./area/4; 
    diag1 = -w12./l2.^2 - w31./l3.^2;
    diag2 = -w12./l1.^2 - w23./l3.^2;
    diag3 = -w23./l2.^2 - w31./l1.^2;

    II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
    JJ = [f3; f3; f1; f1; f2; f2; f1; f2; f3];
    V = [w12./(l2.^2); w12./(l1.^2); w23./(l3.^2); w23./(l2.^2); w31./(l1.^2); w31./(l3.^2); diag1; diag2; diag3];
    M = sparse(II,JJ,V,nv,nv);
    
elseif strcmp(weight,'tutte')
    % The Tutte weight (M_ij = 1)
    nf = length(f);

    I = reshape(f',nf*3,1);
    J = reshape(f(:,[2 3 1])',nf*3,1);
    V = ones(nf*3,1)/2;
    W = sparse([I;J],[J;I],[V;V]);
    M = W + sparse(1:nv,1:nv,-diag(W)-(sum(W)'), nv, nv);

elseif strcmp(weight,'squared_inverse_distance')
    % The squared inverse distance weight (M_ij = 1/|v_i - v_j|^2)
    f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

    l1 = sqrt(sum((v(f2,:) - v(f3,:)).^2,2));
    l2 = sqrt(sum((v(f3,:) - v(f1,:)).^2,2));
    l3 = sqrt(sum((v(f1,:) - v(f2,:)).^2,2));

    w12 = 1./l3.^2/2;
    w23 = 1./l1.^2/2;
    w31 = 1./l2.^2/2;
    diag1 = -w12-w31; diag2 = -w12-w23; diag3 = -w31-w23;

    II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
    JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
    V = [w12; w12; w23; w23; w31; w31; diag1; diag2; diag3];
    M = sparse(II,JJ,V,nv,nv);
  
elseif strcmp(weight,'inverse_distance')
    % The inverse distance weight (M_ij = 1/|v_i - v_j|)
    f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

    l1 = sqrt(sum((v(f2,:) - v(f3,:)).^2,2));
    l2 = sqrt(sum((v(f3,:) - v(f1,:)).^2,2));
    l3 = sqrt(sum((v(f1,:) - v(f2,:)).^2,2));

    w12 = 1./l3/2;
    w23 = 1./l1/2;
    w31 = 1./l2/2;
    diag1 = -w12-w31; diag2 = -w12-w23; diag3 = -w31-w23;

    II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
    JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
    V = [w12; w12; w23; w23; w31; w31; diag1; diag2; diag3];
    M = sparse(II,JJ,V,nv,nv);  
else error('unknown weight!');
end

if normalization == 1
    % normalize
    row_max = abs(diag(M));
    [Mrow,Mcol,Mval] = find(M);
    M = sparse(Mrow,Mcol,Mval./row_max(Mrow),nv,nv);
end