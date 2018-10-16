function [residual, g1, g2, g3] = SmetsBenchmarkEstimation_static(y, x, params)
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

residual = zeros( 28, 1);

%
% Model equations
%

T38 = 1/(1+params(4)*params(2));
T41 = params(2)^2;
T66 = params(28)/params(2);
T71 = (1-T66)/(params(34)*(1+T66));
T114 = 1/(1+params(4)*params(2)*params(32));
T128 = (1-params(33))*(1-params(4)*params(2)*params(33))/params(33)/(1+(params(29)-1)*params(17));
T136 = params(4)*params(2)/(1+params(4)*params(2));
T163 = (1-params(31))*(1-params(4)*params(2)*params(31))/((1+params(4)*params(2))*params(31))*1/(1+(params(20)-1)*params(16));
lhs =y(8);
rhs =params(24)*y(10)+(1-params(24))*y(18)-y(21);
residual(1)= lhs-rhs;
lhs =y(9);
rhs =y(10)*(1-params(26))/params(26);
residual(2)= lhs-rhs;
lhs =y(10);
rhs =y(18)+y(16)-y(11);
residual(3)= lhs-rhs;
lhs =y(11);
rhs =y(9)+y(28);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =T38*(y(14)+y(14)*params(4)*params(2)+1/(T41*params(27))*y(12))+y(24);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =y(12)*(1-params(19))/(1-params(19)+params(5))+y(10)*params(5)/(1-params(19)+params(5))-y(19)+y(17)+1/T71*y(22);
residual(6)= lhs-rhs;
lhs =y(13);
rhs =y(22)+y(13)*T66/(1+T66)+y(13)*1/(1+T66)-T71*(y(19)-y(17));
residual(7)= lhs-rhs;
lhs =y(15);
rhs =y(13)*params(12)+y(14)*params(11)+y(23)+y(9)*params(13);
residual(8)= lhs-rhs;
lhs =y(15);
rhs =params(29)*(y(21)+params(24)*y(11)+(1-params(24))*y(16));
residual(9)= lhs-rhs;
lhs =y(17);
rhs =T114*(params(4)*params(2)*y(17)+y(17)*params(32)+y(8)*T128)+y(26);
residual(10)= lhs-rhs;
lhs =y(18);
rhs =y(18)*T38+y(18)*T136+y(17)*params(30)/(1+params(4)*params(2))-y(17)*(1+params(4)*params(2)*params(30))/(1+params(4)*params(2))+y(17)*T136+T163*(y(16)*params(35)+y(13)*1/(1-T66)-y(13)*T66/(1-T66)-y(18))+y(27);
residual(11)= lhs-rhs;
lhs =y(19);
rhs =y(17)*params(36)*(1-params(39))+(1-params(39))*params(38)*(y(15)-y(21)*params(29))+y(19)*params(39)+y(25);
residual(12)= lhs-rhs;
lhs =y(21);
rhs =y(21)*params(40)+x(1);
residual(13)= lhs-rhs;
lhs =y(23);
rhs =y(23)*params(47)+x(3)+x(1)*params(46);
residual(14)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(44)+x(4);
residual(15)= lhs-rhs;
lhs =y(22);
rhs =y(22)*params(41)+x(2);
residual(16)= lhs-rhs;
lhs =y(26);
rhs =y(26)*params(42)+x(6)-x(6)*params(49);
residual(17)= lhs-rhs;
lhs =y(27);
rhs =y(27)*params(43)+x(7)-x(7)*params(48);
residual(18)= lhs-rhs;
lhs =y(25);
rhs =y(25)*params(45)+x(5);
residual(19)= lhs-rhs;
lhs =y(28);
rhs =y(28)*(1-params(7))+y(14)*params(7)+y(24)*params(27)*T41*params(7);
residual(20)= lhs-rhs;
lhs =y(20);
rhs =y(15)-y(21)*params(29);
residual(21)= lhs-rhs;
lhs =y(4);
rhs =params(25);
residual(22)= lhs-rhs;
lhs =y(5);
rhs =params(25);
residual(23)= lhs-rhs;
lhs =y(6);
rhs =params(25);
residual(24)= lhs-rhs;
lhs =y(7);
rhs =params(25);
residual(25)= lhs-rhs;
lhs =y(3);
rhs =y(17)+params(22);
residual(26)= lhs-rhs;
lhs =y(2);
rhs =y(19)+params(15);
residual(27)= lhs-rhs;
lhs =y(1);
rhs =y(16)+params(21);
residual(28)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(28, 28);

  %
  % Jacobian matrix
  %

  g1(1,8)=1;
  g1(1,10)=(-params(24));
  g1(1,18)=(-(1-params(24)));
  g1(1,21)=1;
  g1(2,9)=1;
  g1(2,10)=(-((1-params(26))/params(26)));
  g1(3,10)=1;
  g1(3,11)=1;
  g1(3,16)=(-1);
  g1(3,18)=(-1);
  g1(4,9)=(-1);
  g1(4,11)=1;
  g1(4,28)=(-1);
  g1(5,12)=(-(T38*1/(T41*params(27))));
  g1(5,14)=1-(1+params(4)*params(2))*T38;
  g1(5,24)=(-1);
  g1(6,10)=(-(params(5)/(1-params(19)+params(5))));
  g1(6,12)=1-(1-params(19))/(1-params(19)+params(5));
  g1(6,17)=(-1);
  g1(6,19)=1;
  g1(6,22)=(-(1/T71));
  g1(7,13)=1-(T66/(1+T66)+1/(1+T66));
  g1(7,17)=(-T71);
  g1(7,19)=T71;
  g1(7,22)=(-1);
  g1(8,9)=(-params(13));
  g1(8,13)=(-params(12));
  g1(8,14)=(-params(11));
  g1(8,15)=1;
  g1(8,23)=(-1);
  g1(9,11)=(-(params(24)*params(29)));
  g1(9,15)=1;
  g1(9,16)=(-((1-params(24))*params(29)));
  g1(9,21)=(-params(29));
  g1(10,8)=(-(T114*T128));
  g1(10,17)=1-T114*(params(4)*params(2)+params(32));
  g1(10,26)=(-1);
  g1(11,13)=(-(T163*(1/(1-T66)-T66/(1-T66))));
  g1(11,16)=(-(T163*params(35)));
  g1(11,17)=(-(T136+params(30)/(1+params(4)*params(2))-(1+params(4)*params(2)*params(30))/(1+params(4)*params(2))));
  g1(11,18)=1-(T38+T136-T163);
  g1(11,27)=(-1);
  g1(12,15)=(-((1-params(39))*params(38)));
  g1(12,17)=(-(params(36)*(1-params(39))));
  g1(12,19)=1-params(39);
  g1(12,21)=(-((1-params(39))*params(38)*(-params(29))));
  g1(12,25)=(-1);
  g1(13,21)=1-params(40);
  g1(14,23)=1-params(47);
  g1(15,24)=1-params(44);
  g1(16,22)=1-params(41);
  g1(17,26)=1-params(42);
  g1(18,27)=1-params(43);
  g1(19,25)=1-params(45);
  g1(20,14)=(-params(7));
  g1(20,24)=(-(params(27)*T41*params(7)));
  g1(20,28)=1-(1-params(7));
  g1(21,15)=(-1);
  g1(21,20)=1;
  g1(21,21)=params(29);
  g1(22,4)=1;
  g1(23,5)=1;
  g1(24,6)=1;
  g1(25,7)=1;
  g1(26,3)=1;
  g1(26,17)=(-1);
  g1(27,2)=1;
  g1(27,19)=(-1);
  g1(28,1)=1;
  g1(28,16)=(-1);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],28,784);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],28,21952);
end
end
end
end
