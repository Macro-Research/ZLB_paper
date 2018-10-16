function [S,P]=kalman_init(gamma1,gamma2,gamma3,sigma);
vec_Sigma=reshape(sigma,[length(sigma)^2,1]);
vec_P= (eye(length(gamma1)^2)-kron(gamma1,gamma1))^(-1)*kron(gamma3,gamma3)*vec_Sigma;
size=sqrt(length(vec_P));
P=reshape(vec_P,[size,size]);
S=(eye(length(gamma1))-gamma1)^(-1)*gamma2;
end