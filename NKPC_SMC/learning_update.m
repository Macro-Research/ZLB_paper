function [ alpha,beta] = learning_update( X,T)

lengthX = length(X);




alpha= sum(X)/T;

beta_num=0;
beta_denum=0;
 for i=1:T-1
      beta_num= beta_num+(X(i)-alpha)*(X(i+1)-alpha);
  end

for i=1:T
  beta_denum = beta_denum + (X(i)-alpha)^2 ; 
end
  
  beta = beta_num/beta_denum;

    
end

