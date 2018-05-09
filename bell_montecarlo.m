%% Monte-Carlo simulation of Bell test
% DKS
% 2018-05-08


%% config
n_shots=100;
Nsc0=10;
dk_unc=[0,0,0];        
det_qe = 1;
det_ring=0;
n_background=0;

theta0=0;
n_rot_axis=[0,1,0];     % rotate around Y-axis

verbose=1;


v_src=1/sqrt(2)*[0 1 1 0]';     % anti-corr triplet


%% main
%%% scattering halo
% generate back-to-back correlated source @ 100% detector QE
%   NOTE: k_A/k_B are momenta at modes A/B and NOT spins.
[k_A,k_B]=genCorrSource(Nsc0,dk_unc,1,verbose);



%%% Bloch sphere
% get rotated 2-qubit state
v_rot=kron(R_n(n_rot_axis,theta0),R_n(n_rot_axis,theta0))*v_src;
p_rot=densitymatrix(v_rot);

% measure spin in original basis
idx_basis_meas=measureDM(p_rot,Nsc0);
% NOTE:
%   canonical basis:
%       1: A-up; B-up;
%       2: A-up; B-down;
%       3: A-down; B-up;
%       4: A-down; B-down;


%%% build halo - spin and momentum
k_1=
k_2=




 
% 
% % randomise count ordering
% idx_rand_ord=randperm(size(k_B,1));      % suffices to randomise k2
% k_B=k_B(idx_rand_ord,:);
% 
% 
