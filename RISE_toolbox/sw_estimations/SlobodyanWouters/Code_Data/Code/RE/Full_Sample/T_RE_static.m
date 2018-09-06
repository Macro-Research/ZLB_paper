function [residual, g1, g2, g3] = T_RE_static(y, x, params)
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

residual = zeros( 36, 1);

%
% Model equations
%

T23 = 1/(params(10)/(1-params(10)));
T39 = 1/(1+params(44)*params(42));
T42 = params(42)^2;
T59 = params(15)/params(42);
T64 = (1-T59)/(params(14)*(1+T59));
T116 = 1/(1+params(44)*params(42)*params(21));
T130 = (1-params(22))*(1-params(44)*params(42)*params(22))/params(22)/(1+(params(18)-1)*params(3));
T138 = params(44)*params(42)/(1+params(44)*params(42));
T165 = (1-params(20))*(1-params(44)*params(42)*params(20))/((1+params(44)*params(42))*params(20))*1/(1+(params(24)-1)*params(1));
lhs =y(10);
rhs =params(9)*y(12)+(1-params(9))*y(20)-y(22);
residual(1)= lhs-rhs;
lhs =y(11);
rhs =y(12)*T23;
residual(2)= lhs-rhs;
lhs =y(12);
rhs =y(20)+y(18)-y(13);
residual(3)= lhs-rhs;
lhs =y(13);
rhs =y(11)+y(29);
residual(4)= lhs-rhs;
lhs =y(16);
rhs =T39*(y(16)+y(16)*params(44)*params(42)+1/(T42*params(12))*y(14))+y(25);
residual(5)= lhs-rhs;
lhs =y(14);
rhs =(-y(21))+y(19)+y(23)*1/T64+y(12)*params(47)/(params(47)+1-params(13))+y(14)*(1-params(13))/(params(47)+1-params(13));
residual(6)= lhs-rhs;
lhs =y(15);
rhs =y(23)+y(15)*T59/(1+T59)+y(15)*1/(1+T59)-T64*(y(21)-y(19));
residual(7)= lhs-rhs;
lhs =y(17);
rhs =y(15)*params(54)+y(16)*params(53)+y(24)+y(11)*params(55);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =params(18)*(y(22)+params(9)*y(13)+(1-params(9))*y(18));
residual(9)= lhs-rhs;
lhs =y(19);
rhs =T116*(params(44)*params(42)*y(19)+y(19)*params(21)+y(10)*T130)+y(27);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =y(20)*T39+y(20)*T138+y(19)*params(19)/(1+params(44)*params(42))-y(19)*(1+params(44)*params(42)*params(19))/(1+params(44)*params(42))+y(19)*T138+T165*(y(18)*params(23)+y(15)*1/(1-T59)-y(15)*T59/(1-T59)-y(20))+y(28);
residual(11)= lhs-rhs;
lhs =y(21);
rhs =y(19)*params(26)*(1-params(29))+(1-params(29))*params(28)*(y(17)-y(22)*params(18))+y(21)*params(29)+y(26);
residual(12)= lhs-rhs;
lhs =y(29);
rhs =y(29)*(1-params(49))+y(16)*params(49)+y(25)*params(12)*T42*params(49);
residual(13)= lhs-rhs;
lhs =y(22);
rhs =y(22)*params(30)+x(1);
residual(14)= lhs-rhs;
lhs =y(23);
rhs =y(23)*params(32)+x(2);
residual(15)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(33)+x(3)+x(1)*params(2);
residual(16)= lhs-rhs;
lhs =y(25);
rhs =y(25)*params(35)+x(4);
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(26)*params(36)+x(5);
residual(18)= lhs-rhs;
lhs =y(27);
rhs =y(27)*params(37)+y(9)-y(9)*params(8);
residual(19)= lhs-rhs;
lhs =y(9);
rhs =x(6);
residual(20)= lhs-rhs;
lhs =y(28);
rhs =y(28)*params(38)+y(8)-y(8)*params(7);
residual(21)= lhs-rhs;
lhs =y(8);
rhs =x(7);
residual(22)= lhs-rhs;
lhs =y(4);
rhs =params(39);
residual(23)= lhs-rhs;
lhs =y(5);
rhs =params(39);
residual(24)= lhs-rhs;
lhs =y(6);
rhs =params(39);
residual(25)= lhs-rhs;
lhs =y(7);
rhs =params(39);
residual(26)= lhs-rhs;
lhs =y(3);
rhs =y(19)+params(5);
residual(27)= lhs-rhs;
lhs =y(2);
rhs =y(21)+params(40);
residual(28)= lhs-rhs;
lhs =y(1);
rhs =y(18)+params(4);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =y(15);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =y(16);
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(18);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =y(19);
residual(33)= lhs-rhs;
lhs =y(34);
rhs =y(14);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =y(12);
residual(35)= lhs-rhs;
lhs =y(36);
rhs =y(20);
residual(36)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(36, 36);

  %
  % Jacobian matrix
  %

  g1(1,10)=1;
  g1(1,12)=(-params(9));
  g1(1,20)=(-(1-params(9)));
  g1(1,22)=1;
  g1(2,11)=1;
  g1(2,12)=(-T23);
  g1(3,12)=1;
  g1(3,13)=1;
  g1(3,18)=(-1);
  g1(3,20)=(-1);
  g1(4,11)=(-1);
  g1(4,13)=1;
  g1(4,29)=(-1);
  g1(5,14)=(-(T39*1/(T42*params(12))));
  g1(5,16)=1-(1+params(44)*params(42))*T39;
  g1(5,25)=(-1);
  g1(6,12)=(-(params(47)/(params(47)+1-params(13))));
  g1(6,14)=1-(1-params(13))/(params(47)+1-params(13));
  g1(6,19)=(-1);
  g1(6,21)=1;
  g1(6,23)=(-(1/T64));
  g1(7,15)=1-(T59/(1+T59)+1/(1+T59));
  g1(7,19)=(-T64);
  g1(7,21)=T64;
  g1(7,23)=(-1);
  g1(8,11)=(-params(55));
  g1(8,15)=(-params(54));
  g1(8,16)=(-params(53));
  g1(8,17)=1;
  g1(8,24)=(-1);
  g1(9,13)=(-(params(9)*params(18)));
  g1(9,17)=1;
  g1(9,18)=(-((1-params(9))*params(18)));
  g1(9,22)=(-params(18));
  g1(10,10)=(-(T116*T130));
  g1(10,19)=1-T116*(params(44)*params(42)+params(21));
  g1(10,27)=(-1);
  g1(11,15)=(-(T165*(1/(1-T59)-T59/(1-T59))));
  g1(11,18)=(-(T165*params(23)));
  g1(11,19)=(-(T138+params(19)/(1+params(44)*params(42))-(1+params(44)*params(42)*params(19))/(1+params(44)*params(42))));
  g1(11,20)=1-(T39+T138-T165);
  g1(11,28)=(-1);
  g1(12,17)=(-((1-params(29))*params(28)));
  g1(12,19)=(-(params(26)*(1-params(29))));
  g1(12,21)=1-params(29);
  g1(12,22)=(-((1-params(29))*params(28)*(-params(18))));
  g1(12,26)=(-1);
  g1(13,16)=(-params(49));
  g1(13,25)=(-(params(12)*T42*params(49)));
  g1(13,29)=1-(1-params(49));
  g1(14,22)=1-params(30);
  g1(15,23)=1-params(32);
  g1(16,24)=1-params(33);
  g1(17,25)=1-params(35);
  g1(18,26)=1-params(36);
  g1(19,9)=(-(1-params(8)));
  g1(19,27)=1-params(37);
  g1(20,9)=1;
  g1(21,8)=(-(1-params(7)));
  g1(21,28)=1-params(38);
  g1(22,8)=1;
  g1(23,4)=1;
  g1(24,5)=1;
  g1(25,6)=1;
  g1(26,7)=1;
  g1(27,3)=1;
  g1(27,19)=(-1);
  g1(28,2)=1;
  g1(28,21)=(-1);
  g1(29,1)=1;
  g1(29,18)=(-1);
  g1(30,15)=(-1);
  g1(30,30)=1;
  g1(31,16)=(-1);
  g1(31,31)=1;
  g1(32,18)=(-1);
  g1(32,32)=1;
  g1(33,19)=(-1);
  g1(33,33)=1;
  g1(34,14)=(-1);
  g1(34,34)=1;
  g1(35,12)=(-1);
  g1(35,35)=1;
  g1(36,20)=(-1);
  g1(36,36)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],36,1296);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],36,46656);
end
end
end
end
