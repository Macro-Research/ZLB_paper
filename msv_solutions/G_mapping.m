function[GG] = G_mapping(beta,gamma1,gamma2,gamma3,varCovar_vec,numVar)


%beta=diag(diag(beta));
MM=gamma1+gamma2*beta^2;

vec0=(eye(numVar^2)- kron(MM,MM))^(-1)*...
    (kron(gamma3,gamma3))*varCovar_vec;

vec1= (kron(eye(numVar),gamma1)+kron(eye(numVar),gamma2*beta^2))*vec0;

GG = reshape(vec1./vec0,[numVar,numVar]);

 GG(find(isnan(GG)))=0; 
 GG_aux=GG;
 GG=zeros(5,5);
 GG(1:3,1:3)=GG_aux(1:3,1:3);
%  GG(3:5,:)=zeros(3,numVar);
 
 %GG=diag(diag(GG));

%fixedPoint= G-beta;



end




