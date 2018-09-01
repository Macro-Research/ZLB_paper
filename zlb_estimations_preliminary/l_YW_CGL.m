function[alpha_tt,beta_tt]=l_YW_CGL(X,alphaOld,gain,numVar)

tt=size(X,2);
 alpha_tt=alphaOld+gain*(X(:,tt)-alphaOld);
    
 Xhat=X(:,1:tt)-alpha_tt*ones(1,tt);
    
    
    Zhat=zeros(numVar,numVar,tt);
    Vhat=zeros(numVar,numVar,tt);
    
    Zhat(:,:,2)=gain*Xhat(:,2)*Xhat(:,1)';
    Vhat(:,:,2)=gain*Xhat(:,2)*Xhat(:,2)'+gain*(1-gain)*Xhat(:,1)*Xhat(:,1)';
    
     for jj=3:tt
         Zhat(:,:,jj)=Zhat(:,:,jj-1)+gain*(Xhat(:,jj)*Xhat(:,jj-1)'-Zhat(:,:,jj-1));
         Vhat(:,:,jj)=Vhat(:,:,jj-1)+gain*(Xhat(:,jj)*Xhat(:,jj)'-Vhat(:,:,jj-1));
     end
      beta_tt=sqrt(1-gain)*Zhat(:,:,tt)*Vhat(:,:,tt)^(-1);
%      beta_tt(beta_tt>2)=0;
%      beta_tt(beta_tt<-2)=0;

end