function U = R_n(n_axis,alpha)
    % Qubit Bloch sphere rotation about an arbitrary axis
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