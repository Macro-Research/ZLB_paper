function [lik] =sigmaPointFilter_NKPC(param,dataset)
param1=[param(1),param(2),param(3),param(4),param(5),param(6),param(7),...
    param(8),param(9),param(10),param(11),param(12),param(13)];
param2=[param(1),param(2),param(3),param(4),param(5),param(14),param(15),...
    param(8),param(9),param(16),param(17),param(18),param(19)];
p_LN=param(20);p_NL=param(21); 
% p_LN=0.1;p_NL=0.05;

% param1=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
% param2=[0 0 0 0.01 3 0 0 0.5 0.5 0 0.07 0.03 0.03 ];
% load('simulated_dataset.mat');

size_dataset=size(dataset);l=size_dataset(2);
sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime2(param2);

numVar=5;

A1_inv=A1^(-1);A2_inv=A2^(-1);


Q=[1-p_LN,p_LN;p_NL,1-p_NL];
alpha1=zeros(numVar,1);
beta1=0.5*eye(numVar);

T= size(dataset,1);
lik=zeros(T,1);


Sigma1=diag([sigma_y1^2;sigma_pinf1^2;sigma_r1^2]);
Sigma2=diag([sigma_y2^2;sigma_pinf2^2;sigma_r2^2]);


gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A1_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A1_inv*D2;
S0=zeros(numVar,1);P0=0.1*eye(numVar);

H=0.0*diag(var(dataset));
s_fore1=zeros(numVar,1);
s_fore2=zeros(numVar,1);
p_fore1=eye(numVar);
p_fore2=eye(numVar);
prob_regime=0.5*ones(2,T);
prob_regime(:,1)=[0.5 0.5];

prob_update=0.5*ones(2,T);

s_up1=zeros(T,numVar);
s_up2=zeros(T,numVar);

numSigmaPoints=2*numVar+1;
sigmaWeights=ones(numSigmaPoints,2)/numSigmaPoints;
sHat_sigma=zeros(2,numVar,numSigmaPoints);
sHat=zeros(2,numVar);
pHat(1,:,:)=zeros(numVar);pHat(2,:,:)=zeros(numVar);
pp=nan(2,2);
for tt=1:1:T
    %forecasting y_t & forecast errors
    yy = dataset(tt,:);
    v1  = yy'-E1 - F1*s_fore1;
    v2  = yy'-E2 - F2*s_fore2;
    Fe1  = F1*p_fore1*F1' + H;
    Fe2  = F2*p_fore2*F2' + H;
    
    %marginal likelihoods
    ml1=-0.5*l*log(2*pi)-0.5*log(det(Sigma1))-0.5*v1'*((Sigma1)\v1);
    ml2=-0.5*l*log(2*pi)-0.5*log(det(Sigma2))-0.5*v2'*((Sigma2)\v2);
    
    lik(tt)=exp(ml1)*prob_regime(1,tt)+exp(ml2)*prob_regime(2,tt);%prob_regime=p_t|{t-1}
    
    %updating of states
    prob_update(1,tt)=(exp(ml1)*prob_regime(1,tt))/(lik(tt));%prob_update=p_t|t
    prob_update(2,tt)=(exp(ml2)*prob_regime(2,tt))/(lik(tt));%prob_update=p_t|t
    
%     p_XY1=p_fore1*F1';
%     p_XY2=p_fore2*F2';
p_XY1=F1*p_fore1;p_XY2=F2*p_fore2;
kGain1=p_fore1*F1'*(F1*p_fore1*F1')^(-1);
kGain2=p_fore2*F2'*(F2*p_fore1*F2')^(-1);
%     kGain1=p_XY1*Fe1^(-1);
%     kGain2=p_XY2*Fe2^(-1);
    
    s_up1(tt,:)=s_fore1+kGain1*v1;
    s_up2(tt,:)=s_fore2+kGain2*v2;
    p_up1(tt,:,:)=p_fore1-kGain1*(F1*p_fore1);
    p_up2(tt,:,:)=p_fore2-kGain2*(F2*p_fore2);
%     p_up1=p_fore1-kGain1*p_XY1';
%     p_up2=p_fore2-kGain2*p_XY2';
    
    %forecasting states
     prob_regime(:,tt+1)=Q'*[prob_update(1,tt);prob_update(2,tt)];
    
     
     for rt1=[1 2]
         for rt=[1 2]
     pp(rt,rt1)=(Q(rt,rt1)*prob_update(rt,tt))/prob_regime(rt1,tt+1);%pp=p(rt,rt+1)
         end
% %          collapsing:


%      s_fore1=gamma1_1*s_up1(tt,:)';
%      s_fore2=gamma1_2*s_up2(tt,:)';
%      p_fore1=gamma1_1*p_up1*gamma1_1'+gamma3_1*Sigma1*gamma3_1';
%      p_fore2=gamma1_2*p_up2*gamma1_2'+gamma3_2*Sigma2*gamma3_2';
%      
      sHat(rt1,:)=pp(1,rt1)*s_up1(tt,:)'+pp(2,rt1)*s_up2(tt,:)';
      pHat(rt1,:,:)=pp(1,rt1)*(gamma1_1*squeeze(p_up1(tt,:,:))*gamma1_1'+gamma3_1*Sigma1*gamma3_1')+...
      +pp(2,rt1)*(gamma1_2*squeeze(p_up2(tt,:,:))*gamma1_2'+gamma3_2*Sigma2*gamma3_2');
      pHat(rt1,:,:)=nearestSPD(squeeze(pHat(rt1,:,:)));
     end
        
     for rt1=1:2
% 
% %sigma point generation
      if rt1==1
          tempP=(squeeze(pHat(rt1,:,:)));
          
      [temp1,temp2]=sigmaPointGenerator(sHat(rt1,:)',tempP...
          ,gamma1_1,gamma2_1,gamma3_1,Sigma1);
      sHat_sigma(rt1,:,:)=temp1';sigmaWeights(:,rt1)=temp2';
      elseif rt1==2
          tempP=(squeeze(pHat(rt1,:,:)));
          [temp1,temp2]=sigmaPointGenerator(sHat(rt1,:)',tempP...
          ,gamma1_2,gamma2_2,gamma3_2,Sigma2);
      sHat_sigma(rt1,:,:)=temp1';sigmaWeights(:,rt1)=temp2';
      end
     end

         
   
     
%      s_fore1=sHat(1,:)';
%      s_fore2=sHat(2,:)';
     p_fore1=squeeze(pHat(1,:,:));
     p_fore2=squeeze(pHat(2,:,:));
     
    s_fore1=zeros(5,1);s_fore2=zeros(5,1);
%     p_fore1=zeros(5,5);p_fore2=zeros(5,5);
    
     for sss=1:numSigmaPoints
         sHat_aux1=reshape(sHat_sigma(1,:,sss),[numVar 1]);
         sHat_aux2=reshape(sHat_sigma(2,:,sss),[numVar 1]);
     s_fore1=s_fore1+sigmaWeights(sss,1)*(gamma1_1*sHat_aux1);
     s_fore2=s_fore2+sigmaWeights(sss,2)*(gamma1_2*sHat_aux2);
     end
%      
%       for sss=1:numSigmaPoints
%          sHat_aux1=reshape(sHat_sigma(1,:,sss),[numVar 1]);
%          sHat_aux2=reshape(sHat_sigma(2,:,sss),[numVar 1]);
% %           p_fore1=p_fore1+sigmaWeights(sss,1)*(((sHat_aux1-s_fore1)*(sHat_aux1-s_fore1).'));
% %           p_fore2=p_fore2+sigmaWeights(sss,2)*(((sHat_aux2-s_fore2)*(sHat_aux2-s_fore2).'));
% %          p_fore1=nearestSPD(p_fore1);p_fore2=nearestSPD(p_fore2);
% %    p_fore1=p_fore1+sigmaWeights(sss,1)* ( (sHat_aux1-s_fore1)*(sHat_aux1-s_fore1)');   
% %    p_fore2=p_fore2+sigmaWeights(sss,2)* ( (sHat_aux2-s_fore2)*(sHat_aux2-s_fore2)');   
% 
%       end
     p_fore1=nearestSPD(p_fore1);p_fore2=nearestSPD(p_fore2);
    
     
    
end
lik=-(sum(log(lik(5:end))));