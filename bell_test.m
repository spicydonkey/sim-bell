%% testing ground for simulating Bell experiment
% DKS
% 2018-05-03


%% config
%%% useful structures
% Pauli matrices
pauliS={
    [0  1;     
    1   0];
    [0  -1i;     
    1i   0];
    [1  0;     
    0   -1];
        };


%%% graphics
cc=distinguishable_colors(4);
lstyle={'-','--','-.',':'};




%% main
%%% zero rotation --> identity
% R_x(0)
% R_y(0)
% R_z(0)


%%% Kronecker product
theta0=pi/2;

theta_A=theta0;
theta_B=theta0;

Uyy=kron(R_y(theta_A),R_y(theta_B));


%%% Bell state 
%   represented in canonical basis: {|+,+>, |+,->, |-,+>,|-,->}
v_bell={1/sqrt(2)*[0 1 1 0]';
        1/sqrt(2)*[0 1 -1 0]';
        1/sqrt(2)*[1 0 0 1]';
        1/sqrt(2)*[1 0 0 -1]'};
bell_namestr={'Psi+','Psi-','Phi+','Phi-'};


v_bell_psi_plus=v_bell{1};
v_bell_psi_min=v_bell{2};
v_bell_phi_plus=v_bell{3};
v_bell_phi_min=v_bell{4};




%% Example 1: basics
%   take some state
% V=v_bell_psi_min
V=v_bell_psi_plus;
% V=v_bell_phi_plus
% V=v_bell_phi_min


% try a rotation transform
% V=Uyy*V
V=Ryy(0,0)*V;


% evaluate density matrix
% pm=V*(V')
pm=densitymatrix(V);


% get measurement probabilities of spin in Z-basis
p=diag(pm);


%%% spin-correlations
corr_spin=[1 -1 -1 1]';      % correlation coefficients for the basis states
Ocorr_spin=diag([1,-1,-1,1]);     % spin-correlation operator


% measure spin-correlation
% NOTE: is this strictly correct for any operator in matrix rep?
E_avg=trace(Ocorr_spin*pm);      % NOTE: works ONLY because Ocorr_spin is diagnonal in this basis



%% Example 2: basis(theta)-dependency of pair spin-correlation
%   here we assume rotation about X/Y-axis then measurement in Z-basis is
%   equivalent to a measurement of the original state in a rotated basis

th0=linspace(-pi,pi,1e3);
E0=cell(4,1);

% evaluate for each bell state
for ii=1:4
    tv=v_bell{ii};      % get this bell state
    
    tE0=NaN(size(th0));
    for jj=1:length(th0)
        tv_rot=Ryy(th0(jj),th0(jj))*tv;
%         tv_rot=Rxx(th0(jj),th0(jj))*tv;
        
        
%         % Method1
%         tpm=densitymatrix(tv_rot);
%         tE0(jj)=trace(Ocorr_spin*tpm);
        
        % Method 2
        tE0(jj)=corr_spin'*(tv_rot.*conj(tv_rot));
        
    end
    E0{ii}=tE0;
end

%%% vis
h=figure('Name','bell_basis_dep');
p=[];
for ii=1:4
    hold on;
    p(ii)=plot(th0/pi,E0{ii},'DisplayName',bell_namestr{ii},...
        'LineWidth',1.5,'Color',cc(ii,:),'LineStyle',lstyle{ii});
end

% annotate
title('Bell states: basis dependency of spin-correlations');
xlabel('$\theta/\pi$');
ylabel('$E$');

legend(p);
ylim([-1.1,1.1]);
box on;



%% Example 3: basis(theta)-dependency of pair spin-correlation: arbitrary axis

% define rotation axis as 3D cartesian vector
n_rot_axis=[1 0 10];


th0=linspace(-pi,pi,1e3);
E0=cell(4,1);

% evaluate for each bell state
for ii=1:4
    tv=v_bell{ii};      % get this bell state
    
    tE0=NaN(size(th0));
    for jj=1:length(th0)
        tv_rot=kron(R_n(n_rot_axis,th0(jj)),R_n(n_rot_axis,th0(jj)))*tv;
        
%         % Method1
%         tpm=densitymatrix(tv_rot);
%         tE0(jj)=trace(Ocorr_spin*tpm);
        
        % Method 2
        tE0(jj)=corr_spin'*(tv_rot.*conj(tv_rot));
        
    end
    E0{ii}=tE0;
end

%%% vis
h=figure('Name','bell_global_arbaxis');
p=[];
for ii=1:4
    hold on;
    p(ii)=plot(th0/pi,E0{ii},'DisplayName',bell_namestr{ii},...
        'LineWidth',1.5,'Color',cc(ii,:),'LineStyle',lstyle{ii});
end

% annotate
titlestr=sprintf('Global rotation: n=[%0.2g,%0.2g,%0.2g]',n_rot_axis);
title(titlestr);
xlabel('$\theta/\pi$');
ylabel('$E$');

legend(p);
ylim([-1.1,1.1]);
box on;



%% utilities
%%% Rotations on the Bloch sphere
%   http://www.vcpc.univie.ac.at/~ian/hotlist/qc/talks/bloch-sphere-rotations.pdf

function U = R_x(theta)
    U=[ cos(theta/2)        1i*(-sin(theta/2))    ;
        1i*(-sin(theta/2))  cos(theta/2)        ];
end

function U = R_y(theta)
    U=[ cos(theta/2)    -sin(theta/2)   ;
        sin(theta/2)    cos(theta/2)    ];
end

function U = R_z(theta)
    U=[ exp(-1i*theta/2)    0   ;
        0                   exp(1i*theta/2)];
end

% rotation about an arbitrary axis
function U = R_n(n_axis,alpha)
    %   n_axis: axis of rotation
    %       Form as:    R^3 Cartesian vector (x,y,z) 
    %       	or      R^2 polar vector (theta,phi)
    %   alpha:  rotation angle
    %
    %   http://www.vcpc.univie.ac.at/~ian/hotlist/qc/talks/bloch-sphere-rotations.pdf
    
    % Pauli matrices
    pauliS={
    [0  1;     
    1   0];
    [0  -1i;     
    1i   0];
    [1  0;     
    0   -1];
        };
    
    % get normed rotation axis as a Cartesian 3-vector
    nvec=NaN(1,3);
    if length(n_axis)==3
        nvec=n_axis/vnorm(n_axis);      % normalise the 3-vector
    elseif length(n_axis)==2
        [nvec(1),nvec(1),nvec(3)]=sph2cart(n_axis(1),n_axis(2),1);
    else
        error('input n_axis is in an incorrect format.');
    end

    % spin operator
    sig_n=zeros(2,2);
    for ii=1:3
        sig_n = sig_n + nvec(ii)*pauliS{ii};
    end
    
    % rotation operator
    U=cos(alpha/2)*diag([1,1]) - 1i*sin(alpha/2)*sig_n;
    
end


function U = Rxx(th1,th2)
    U=kron(R_x(th1),R_x(th2));
end

function U = Ryy(th1,th2)
    U=kron(R_y(th1),R_y(th2));
end




%%% Density matrix
function pm = densitymatrix(v)
    pm=v*(v');        % density matrix (NOTE: "'" operator complex conjugates!
end