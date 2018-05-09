function idx_basis_out=measureDM(dm,n)
% Measurement given density matrix
%   dm:     density matrix
%   n:      num of repeats
%   idx_basis_out:  index to measurement collapsed basis state

if ~exist('n','var')
    n=1;
end

p_basis=diag(dm);

P=cumsum(p_basis);
rr=rand(n,1);

idx_basis_out=NaN(n,1);
for ii=1:n
    idx_basis_out(ii)=find(P-rr(ii)>0,1);
end

end