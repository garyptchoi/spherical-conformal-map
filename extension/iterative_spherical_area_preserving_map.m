function map = iterative_spherical_area_preserving_map(v,f)
% A method for computing spherical area-preserving map of a genus-0 closed surface
%
% Remark: 
% This code is used only for demonstrating a possible extension of our 
% spherical conformal (or tutte) map to minimize other distortions. 
% The key idea is to use a simple variation of our proposed linear method 
% to obtain an initial spherical map, and then solve an energy minimization
% problem to obtain a more area-preserving map. 
% Note that the computation is not optimized, and it may take longer time 
% than the other linear methods that we proposed  
% (spherical_conformal_map/spherical_tutte_map/spherical_area_preserving_map).
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% 
% Output:
% map: nv x 3 vertex coordinates of the spherical parameterization
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%% set-up
% one may need to change the parameters below to achieve better convergence in the energy minimization scheme
threshold = 1e-6; % in case the final result is not good enough, set a smaller value
max_num_step = 500; % in case the final result is not good enough, set a larger value
dt = 1e-4; % initial guess of time step for line search

% normalize v
v = v*sqrt(4*pi/sum(face_area(f,v)));

%% Step 1: Initialization using a simple variation of our proposed linear method
S = spherical_area_preserving_map(v,f);
% depending on applications, one may consider using other versions of our proposed linear method
% S = spherical_conformal_map(v,f);
% S = spherical_tutte_map(f);

%% Step 2: Local area distortion energy minimization
% We consider minimizing the Chi energy (see Desbrun et al., Eurographics 
% 2002) on the sphere. The energy is related to the local area distortion 
% of the spherical parameterization.

% One may change the matrix below to other mesh Laplacian matrices for
% minimizing other types of energy
M = mesh_laplacian(v,f,'authalic'); 

E_old = evaluate_energy(S,M);
step = 0;
fprintf('Iteration #%i: Energy = %.8f\n', [step,E_old]);

% We use the steepest descent approach to minimize the energy on the sphere
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');

while step < max_num_step
    % Compute the absolute derivative
	dS = -M*S;
    N = compute_vertex_normal(S,f);
    dS_dot_N = sum(dS.*N,2);
    dS_perp = repmat(dS_dot_N,1,3).*N;
    
    % projection of the derivative onto the sphere
    dS = dS - dS_perp;
   
    % update S with line search
    E = @(k) evaluate_energy(S-k*dS,M);
    dt_opt = fminunc(E,dt,options);
    S = S - dS*dt_opt;

    % normalize S to ensure a spherical shape
    S = [S(:,1)-mean(S(:,1)), S(:,2)-mean(S(:,2)), S(:,3)-mean(S(:,3))];
    S = S./repmat(sqrt(sum(S.^2,2)),1,3);

    step = step + 1;
	E_new = evaluate_energy(S,M);
	fprintf('Iteration #%i: Energy = %.8f\n', [step,E_new]);
    if abs(E_new - E_old) < threshold
		break;
    end
    if E_new - E_old > 1e4
		error('The algorithm diverges. Double check the gradient or the step size.');
    end
	E_old = E_new;
    dt = dt_opt;
end
fprintf('Final Energy = %.8f\n', E_new);
map = S;