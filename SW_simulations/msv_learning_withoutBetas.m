function [theta r largest_eig1 largest_eig2 projectionFlag] =msv_learning3(x,regressor,thetaOld,rOld,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_thr)

largest_eig1=0;largest_eig2=0;
numBackward=length(backward_indices);
lastwarn('Success');
lambda=10e-15;
projectionFlag=0;

 yy=regressor;


 
r=rOld+gain*(yy*yy'-rOld);


try
eig_min=min(abs(eigs(r)));
catch
eig_min=10e-15;
end


    if eig_min>lambda
     theta=thetaOld'+gain*(r\yy)*(x-thetaOld*yy)';
    else
    theta=thetaOld'+gain*(((r+lambda*eye(length(r))))\yy)*(x-thetaOld*yy)';
    end
    
    



end



