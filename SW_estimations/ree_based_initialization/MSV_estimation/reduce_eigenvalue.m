function[newm]=reduce_eigenvalue(m,eig_crit,eps)

    [VV,DD]=eig(m);
    DD=diag(DD);
    ii=find(DD>eig_crit);
    DD(ii)=eig_crit-eps;

newm=real(VV*diag(DD)*pinv(VV));


end