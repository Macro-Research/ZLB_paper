function [residual, g1, g2, g3] = SW_Estimation_New_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 29, 1);

%
% Model equations
%

T23 = 1/(params(13)/(1-params(13)));
T39 = 1/(1+params(36)*params(50));
T42 = params(50)^2;
T59 = params(15)/params(50);
T64 = (1-T59)/(params(21)*(1+T59));
T116 = 1/(1+params(36)*params(50)*params(19));
T130 = (1-params(20))*(1-params(36)*params(50)*params(20))/params(20)/(1+(params(16)-1)*params(2));
T138 = params(36)*params(50)/(1+params(36)*params(50));
T165 = (1-params(18))*(1-params(36)*params(50)*params(18))/((1+params(36)*params(50))*params(18))*1/(1+(params(5)-1)*params(1));
lhs =y(8);
rhs =params(9)*y(10)+(1-params(9))*y(18)-y(20);
residual(1)= lhs-rhs;
lhs =y(9);
rhs =y(10)*T23;
residual(2)= lhs-rhs;
lhs =y(10);
rhs =y(18)+y(16)-y(11);
residual(3)= lhs-rhs;
lhs =y(11);
rhs =y(9)+y(27);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =T39*(y(14)+y(14)*params(36)*params(50)+1/(T42*params(14))*y(12))+y(23);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =(-y(19))+y(17)+y(21)*1/T64+y(10)*params(38)/(params(38)+1-params(4))+y(12)*(1-params(4))/(params(38)+1-params(4));
residual(6)= lhs-rhs;
lhs =y(13);
rhs =y(21)+y(13)*T59/(1+T59)+y(13)*1/(1+T59)-T64*(y(19)-y(17));
residual(7)= lhs-rhs;
lhs =y(15);
rhs =y(13)*params(45)+y(14)*params(44)+y(22)+y(9)*params(46);
residual(8)= lhs-rhs;
lhs =y(15);
rhs =params(16)*(y(20)+params(9)*y(11)+(1-params(9))*y(16));
residual(9)= lhs-rhs;
lhs =y(17);
rhs =T116*(params(36)*params(50)*y(17)+y(17)*params(19)+y(8)*T130)+y(25);
residual(10)= lhs-rhs;
lhs =y(18);
rhs =y(18)*T39+y(18)*T138+y(17)*params(17)/(1+params(36)*params(50))-y(17)*(1+params(36)*params(50)*params(17))/(1+params(36)*params(50))+y(17)*T138+T165*(y(16)*params(22)+y(13)*1/(1-T59)-y(13)*T59/(1-T59)-y(18))+y(26);
residual(11)= lhs-rhs;
lhs =y(19);
rhs =y(17)*params(23)*(1-params(26))+(1-params(26))*params(25)*(y(15)-y(20)*params(16))+params(24)*(y(20)*params(16)+y(15)-y(20)*params(16)-y(15))+y(19)*params(26)+y(24);
residual(12)= lhs-rhs;
lhs =y(20);
rhs =y(20)*params(28)+x(1);
residual(13)= lhs-rhs;
lhs =y(21);
rhs =y(21)*params(29)+x(2);
residual(14)= lhs-rhs;
lhs =y(22);
rhs =y(22)*params(34)+x(3)+x(1)*params(3);
residual(15)= lhs-rhs;
lhs =y(23);
rhs =y(23)*params(32)+x(4);
residual(16)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(33)+x(5);
residual(17)= lhs-rhs;
lhs =y(25);
rhs =y(25)*params(30)+x(6)-x(6)*params(12);
residual(18)= lhs-rhs;
lhs =y(26);
rhs =y(26)*params(31)+x(7)-x(7)*params(11);
residual(19)= lhs-rhs;
lhs =y(27);
rhs =y(27)*(1-params(40))+y(14)*params(40)+y(23)*params(14)*T42*params(40);
residual(20)= lhs-rhs;
lhs =y(4);
rhs =params(10);
residual(21)= lhs-rhs;
lhs =y(5);
rhs =params(10);
residual(22)= lhs-rhs;
lhs =y(6);
rhs =params(10);
residual(23)= lhs-rhs;
lhs =y(7);
rhs =params(10);
residual(24)= lhs-rhs;
lhs =y(3);
rhs =y(17)+params(7);
residual(25)= lhs-rhs;
lhs =y(2);
rhs =y(19)+params(49);
residual(26)= lhs-rhs;
lhs =y(1);
rhs =y(16)+params(6);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =x(6);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =x(7);
residual(29)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(29, 29);

  %
  % Jacobian matrix
  %

  g1(1,8)=1;
  g1(1,10)=(-params(9));
  g1(1,18)=(-(1-params(9)));
  g1(1,20)=1;
  g1(2,9)=1;
  g1(2,10)=(-T23);
  g1(3,10)=1;
  g1(3,11)=1;
  g1(3,16)=(-1);
  g1(3,18)=(-1);
  g1(4,9)=(-1);
  g1(4,11)=1;
  g1(4,27)=(-1);
  g1(5,12)=(-(T39*1/(T42*params(14))));
  g1(5,14)=1-(1+params(36)*params(50))*T39;
  g1(5,23)=(-1);
  g1(6,10)=(-(params(38)/(params(38)+1-params(4))));
  g1(6,12)=1-(1-params(4))/(params(38)+1-params(4));
  g1(6,17)=(-1);
  g1(6,19)=1;
  g1(6,21)=(-(1/T64));
  g1(7,13)=1-(T59/(1+T59)+1/(1+T59));
  g1(7,17)=(-T64);
  g1(7,19)=T64;
  g1(7,21)=(-1);
  g1(8,9)=(-params(46));
  g1(8,13)=(-params(45));
  g1(8,14)=(-params(44));
  g1(8,15)=1;
  g1(8,22)=(-1);
  g1(9,11)=(-(params(9)*params(16)));
  g1(9,15)=1;
  g1(9,16)=(-((1-params(9))*params(16)));
  g1(9,20)=(-params(16));
  g1(10,8)=(-(T116*T130));
  g1(10,17)=1-T116*(params(36)*params(50)+params(19));
  g1(10,25)=(-1);
  g1(11,13)=(-(T165*(1/(1-T59)-T59/(1-T59))));
  g1(11,16)=(-(T165*params(22)));
  g1(11,17)=(-(T138+params(17)/(1+params(36)*params(50))-(1+params(36)*params(50)*params(17))/(1+params(36)*params(50))));
  g1(11,18)=1-(T39+T138-T165);
  g1(11,26)=(-1);
  g1(12,15)=(-((1-params(26))*params(25)));
  g1(12,17)=(-(params(23)*(1-params(26))));
  g1(12,19)=1-params(26);
  g1(12,20)=(-((1-params(26))*params(25)*(-params(16))));
  g1(12,24)=(-1);
  g1(13,20)=1-params(28);
  g1(14,21)=1-params(29);
  g1(15,22)=1-params(34);
  g1(16,23)=1-params(32);
  g1(17,24)=1-params(33);
  g1(18,25)=1-params(30);
  g1(19,26)=1-params(31);
  g1(20,14)=(-params(40));
  g1(20,23)=(-(params(14)*T42*params(40)));
  g1(20,27)=1-(1-params(40));
  g1(21,4)=1;
  g1(22,5)=1;
  g1(23,6)=1;
  g1(24,7)=1;
  g1(25,3)=1;
  g1(25,17)=(-1);
  g1(26,2)=1;
  g1(26,19)=(-1);
  g1(27,1)=1;
  g1(27,16)=(-1);
  g1(28,28)=1;
  g1(29,29)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],29,841);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],29,24389);
end
end
end
end
