function [theta r largest_eig1 largest_eig2 projectionFlag] =msv_learning3(x,regressor,thetaOld,rOld,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2)
msgstr=[];
msgid=[];
numBackward=length(backward_indices);
lastwarn('Success');
lambda=10e-8;
projectionFlag=0;
largest_eig1=0;largest_eig2=0;

 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);

% [VV DD]=eig(r);
% DD=diag(DD);
% DD(DD<lambda)=min(DD(DD>lambda));
% r=VV* diag(DD) / VV;


try
eig_min=min(abs(eigs(r)));
catch
eig_min=10e-15;
end


    if eig_min>lambda
     theta=thetaOld'+gain*(r^(-1)*yy)*(x-thetaOld*yy)';
    else
    theta=thetaOld'+gain*(((r+lambda*eye(length(r))))^(-1)*yy)*(x-thetaOld*yy)';
    end
    
    
beta_tt=zeros(numEndo,numEndo);
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);

try
largest_eig1=abs(eigs(beta_tt,1));
ev1=abs(eigs(gamma1_1_tilde,1));
ev2=abs(eigs(gamma1_2_tilde,1));
largest_eig2=max(ev1,ev2);
catch
    largest_eig1=1.01;
    largest_eig2=1.01;
end

[msgstr, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
r=rOld;
theta=thetaOld';
projectionFlag=1;
elseif largest_eig1>1
    r=rOld;
theta=thetaOld';
projectionFlag=1;

elseif largest_eig2>1
    
        r=rOld;
theta=thetaOld';
projectionFlag=1;

end



    



end



