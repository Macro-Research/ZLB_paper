clear;
clc;close all;
tic
seed=round(1000*rand);
rng(seed);
parameters=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
% parameters=[0 0 0 0.01 3  0  0   0.5 0.5 0   0.7 0.3 0.03 ];

sigma_y = parameters(end-2);
sigma_pinf=parameters(end-1);
sigma_r=parameters(end);

[Atotal Btotal Ctotal Dtotal Etotal Ftotal Gtotal]= NKPC_sysmat(parameters);
N=5000;

    

Xtotal(:,1)=rand(5,1);


alpha1=Xtotal(1,1), beta1=-0.5, r1=0;
alpha2=Xtotal(2,1), beta2=-0.5, r2=0;


learning1(1,:)= [alpha1,beta1];
learning2(1,:)= [alpha2,beta2];



alpha=[alpha1;alpha2];
beta=diag([beta1;beta2]);

betaShocks=zeros(2,1);

alphaTotal=[alpha',0,0,0]';
betaTotal=[beta,zeros(2,3);zeros(3,5)]';
       
eps_y = normrnd(0,sigma_y,[N,1]);
eps_pinf = normrnd(0,sigma_pinf,[N,1]);    
eps_r = normrnd(0,sigma_r,[N,1]);

errors=[eps_y eps_pinf eps_r]' ; 

Atotal_inv = Atotal^(-1);
dist=2;
t=2;


 for t=2:N

disp(t);
      Xtotal(:,t) =  Atotal_inv * ( Btotal*Xtotal(:,t-1)+Ctotal*(alphaTotal+betaTotal^2*(Xtotal(:,t-1)-alphaTotal))+...
      Dtotal*errors(:,t)     ) ;
  
 
%  [alpha1,beta1,r1]=recursive_update(Xtotal(1,:),t,alpha1,beta1,r1);
%  [alpha2,beta2,r2]=recursive_update(Xtotal(2,:),t,alpha2,beta2,r2);
 
%  [alpha1,beta1]=learning_update(Xtotal(1,:),t);
%  [alpha2,beta2]=learning_update(Xtotal(2,:),t);

% [alpha1, beta1]=cgl_learning(Xtotal(1,1:t),alpha1);
% [alpha2, beta2]=cgl_learning(Xtotal(2,1:t),alpha2);
%  
alpha=[alpha1;alpha2];
beta=diag([beta1;beta2]);
 
alphaTotal=[alpha',0,0,0]';
betaTotal=[beta,zeros(2,3);zeros(3,5)]';



learning1(t,:)= [alpha1,beta1];
learning2(t,:)= [alpha2,beta2];


     
 end
 

 figure;
 subplot(2,2,1);
 plot(learning1(:,1),'lineWidth',3);
 legend('Expectation: Output Gap Mean');
 subplot(2,2,2);
 plot(learning1(:,2),'lineWidth',3);
 legend('Expectation: Output Gap Persistence');
 subplot(2,2,3);
 plot(learning2(:,1),'lineWidth',3);
 legend('Expectation: Inflation Mean');
 subplot(2,2,4);
 plot(learning2(:,2),'lineWidth',3);
legend('Expectation: Inflation Persistence');
 
 
   
 

