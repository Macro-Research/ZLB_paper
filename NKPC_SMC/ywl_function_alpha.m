function[theta]=ywl_function_alpha(y,theta)

newObs=y(end);
gain=0.05;
theta=gain*newObs+(1-gain)*theta;


end