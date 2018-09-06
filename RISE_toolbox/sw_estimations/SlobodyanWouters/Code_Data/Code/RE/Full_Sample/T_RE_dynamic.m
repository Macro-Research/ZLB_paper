function [residual, g1, g2, g3] = T_RE_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(36, 1);
T23 = 1/(params(10)/(1-params(10)));
T39 = 1/(1+params(44)*params(42));
T44 = params(42)^2;
T61 = params(15)/params(42);
T66 = (1-T61)/(params(14)*(1+T61));
T94 = (params(14)-1)*params(56)/(params(14)*(1+T61));
T127 = 1/(1+params(44)*params(42)*params(21));
T142 = (1-params(22))*(1-params(44)*params(42)*params(22))/params(22)/(1+(params(18)-1)*params(3));
T151 = params(44)*params(42)/(1+params(44)*params(42));
T179 = (1-params(20))*(1-params(44)*params(42)*params(20))/((1+params(44)*params(42))*params(20))*1/(1+(params(24)-1)*params(1));
lhs =y(29);
rhs =params(9)*y(31)+(1-params(9))*y(39)-y(41);
residual(1)= lhs-rhs;
lhs =y(30);
rhs =y(31)*T23;
residual(2)= lhs-rhs;
lhs =y(31);
rhs =y(39)+y(37)-y(32);
residual(3)= lhs-rhs;
lhs =y(32);
rhs =y(30)+y(19);
residual(4)= lhs-rhs;
lhs =y(35);
rhs =T39*(y(6)+params(44)*params(42)*y(59)+1/(T44*params(12))*y(33))+y(44);
residual(5)= lhs-rhs;
lhs =y(33);
rhs =(-y(40))+y(61)+y(42)*1/T66+params(47)/(params(47)+1-params(13))*y(56)+(1-params(13))/(params(47)+1-params(13))*y(57);
residual(6)= lhs-rhs;
lhs =y(34);
rhs =y(42)+T61/(1+T61)*y(5)+1/(1+T61)*y(58)+T94*(y(37)-y(60))-T66*(y(40)-y(61));
residual(7)= lhs-rhs;
lhs =y(36);
rhs =y(34)*params(54)+y(35)*params(53)+y(43)+y(30)*params(55);
residual(8)= lhs-rhs;
lhs =y(36);
rhs =params(18)*(y(41)+params(9)*y(32)+(1-params(9))*y(37));
residual(9)= lhs-rhs;
lhs =y(38);
rhs =T127*(params(44)*params(42)*y(61)+params(21)*y(9)+y(29)*T142)+y(46);
residual(10)= lhs-rhs;
lhs =y(39);
rhs =T39*y(10)+T151*y(62)+y(9)*params(19)/(1+params(44)*params(42))-y(38)*(1+params(44)*params(42)*params(19))/(1+params(44)*params(42))+y(61)*T151+T179*(y(37)*params(23)+y(34)*1/(1-T61)-y(5)*T61/(1-T61)-y(39))+y(47);
residual(11)= lhs-rhs;
lhs =y(40);
rhs =y(38)*params(26)*(1-params(29))+(1-params(29))*params(28)*(y(36)-y(41)*params(18))+params(27)*(y(36)-y(7)-params(18)*(y(41)-y(12)))+params(29)*y(11)+y(45);
residual(12)= lhs-rhs;
lhs =y(48);
rhs =y(19)*(1-params(49))+y(35)*params(49)+y(44)*params(12)*T44*params(49);
residual(13)= lhs-rhs;
lhs =y(41);
rhs =y(12)*params(30)+x(it_, 1);
residual(14)= lhs-rhs;
lhs =y(42);
rhs =params(32)*y(13)+x(it_, 2);
residual(15)= lhs-rhs;
lhs =y(43);
rhs =params(33)*y(14)+x(it_, 3)+x(it_, 1)*params(2);
residual(16)= lhs-rhs;
lhs =y(44);
rhs =params(35)*y(15)+x(it_, 4);
residual(17)= lhs-rhs;
lhs =y(45);
rhs =params(36)*y(16)+x(it_, 5);
residual(18)= lhs-rhs;
lhs =y(46);
rhs =params(37)*y(17)+y(28)-params(8)*y(2);
residual(19)= lhs-rhs;
lhs =y(28);
rhs =x(it_, 6);
residual(20)= lhs-rhs;
lhs =y(47);
rhs =params(38)*y(18)+y(27)-params(7)*y(1);
residual(21)= lhs-rhs;
lhs =y(27);
rhs =x(it_, 7);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(36)-y(7)+params(39);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =params(39)+y(34)-y(5);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =params(39)+y(35)-y(6);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =params(39)+y(39)-y(10);
residual(26)= lhs-rhs;
lhs =y(22);
rhs =y(38)+params(5);
residual(27)= lhs-rhs;
lhs =y(21);
rhs =y(40)+params(40);
residual(28)= lhs-rhs;
lhs =y(20);
rhs =y(37)+params(4);
residual(29)= lhs-rhs;
lhs =y(49);
rhs =y(5);
residual(30)= lhs-rhs;
lhs =y(50);
rhs =y(6);
residual(31)= lhs-rhs;
lhs =y(51);
rhs =y(8);
residual(32)= lhs-rhs;
lhs =y(52);
rhs =y(9);
residual(33)= lhs-rhs;
lhs =y(53);
rhs =y(4);
residual(34)= lhs-rhs;
lhs =y(54);
rhs =y(3);
residual(35)= lhs-rhs;
lhs =y(55);
rhs =y(10);
residual(36)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(36, 69);

  %
  % Jacobian matrix
  %

  g1(1,29)=1;
  g1(1,31)=(-params(9));
  g1(1,39)=(-(1-params(9)));
  g1(1,41)=1;
  g1(2,30)=1;
  g1(2,31)=(-T23);
  g1(3,31)=1;
  g1(3,32)=1;
  g1(3,37)=(-1);
  g1(3,39)=(-1);
  g1(4,30)=(-1);
  g1(4,32)=1;
  g1(4,19)=(-1);
  g1(5,33)=(-(T39*1/(T44*params(12))));
  g1(5,6)=(-T39);
  g1(5,35)=1;
  g1(5,59)=(-(params(44)*params(42)*T39));
  g1(5,44)=(-1);
  g1(6,56)=(-(params(47)/(params(47)+1-params(13))));
  g1(6,33)=1;
  g1(6,57)=(-((1-params(13))/(params(47)+1-params(13))));
  g1(6,61)=(-1);
  g1(6,40)=1;
  g1(6,42)=(-(1/T66));
  g1(7,5)=(-(T61/(1+T61)));
  g1(7,34)=1;
  g1(7,58)=(-(1/(1+T61)));
  g1(7,37)=(-T94);
  g1(7,60)=T94;
  g1(7,61)=(-T66);
  g1(7,40)=T66;
  g1(7,42)=(-1);
  g1(8,30)=(-params(55));
  g1(8,34)=(-params(54));
  g1(8,35)=(-params(53));
  g1(8,36)=1;
  g1(8,43)=(-1);
  g1(9,32)=(-(params(9)*params(18)));
  g1(9,36)=1;
  g1(9,37)=(-((1-params(9))*params(18)));
  g1(9,41)=(-params(18));
  g1(10,29)=(-(T127*T142));
  g1(10,9)=(-(params(21)*T127));
  g1(10,38)=1;
  g1(10,61)=(-(params(44)*params(42)*T127));
  g1(10,46)=(-1);
  g1(11,5)=(-(T179*(-(T61/(1-T61)))));
  g1(11,34)=(-(T179*1/(1-T61)));
  g1(11,37)=(-(T179*params(23)));
  g1(11,9)=(-(params(19)/(1+params(44)*params(42))));
  g1(11,38)=(1+params(44)*params(42)*params(19))/(1+params(44)*params(42));
  g1(11,61)=(-T151);
  g1(11,10)=(-T39);
  g1(11,39)=1-(-T179);
  g1(11,62)=(-T151);
  g1(11,47)=(-1);
  g1(12,7)=params(27);
  g1(12,36)=(-((1-params(29))*params(28)+params(27)));
  g1(12,38)=(-(params(26)*(1-params(29))));
  g1(12,11)=(-params(29));
  g1(12,40)=1;
  g1(12,12)=(-(params(18)*params(27)));
  g1(12,41)=(-((1-params(29))*params(28)*(-params(18))+params(27)*(-params(18))));
  g1(12,45)=(-1);
  g1(13,35)=(-params(49));
  g1(13,44)=(-(params(12)*T44*params(49)));
  g1(13,19)=(-(1-params(49)));
  g1(13,48)=1;
  g1(14,12)=(-params(30));
  g1(14,41)=1;
  g1(14,63)=(-1);
  g1(15,13)=(-params(32));
  g1(15,42)=1;
  g1(15,64)=(-1);
  g1(16,14)=(-params(33));
  g1(16,43)=1;
  g1(16,63)=(-params(2));
  g1(16,65)=(-1);
  g1(17,15)=(-params(35));
  g1(17,44)=1;
  g1(17,66)=(-1);
  g1(18,16)=(-params(36));
  g1(18,45)=1;
  g1(18,67)=(-1);
  g1(19,2)=params(8);
  g1(19,28)=(-1);
  g1(19,17)=(-params(37));
  g1(19,46)=1;
  g1(20,28)=1;
  g1(20,68)=(-1);
  g1(21,1)=params(7);
  g1(21,27)=(-1);
  g1(21,18)=(-params(38));
  g1(21,47)=1;
  g1(22,27)=1;
  g1(22,69)=(-1);
  g1(23,23)=1;
  g1(23,7)=1;
  g1(23,36)=(-1);
  g1(24,24)=1;
  g1(24,5)=1;
  g1(24,34)=(-1);
  g1(25,25)=1;
  g1(25,6)=1;
  g1(25,35)=(-1);
  g1(26,26)=1;
  g1(26,10)=1;
  g1(26,39)=(-1);
  g1(27,22)=1;
  g1(27,38)=(-1);
  g1(28,21)=1;
  g1(28,40)=(-1);
  g1(29,20)=1;
  g1(29,37)=(-1);
  g1(30,5)=(-1);
  g1(30,49)=1;
  g1(31,6)=(-1);
  g1(31,50)=1;
  g1(32,8)=(-1);
  g1(32,51)=1;
  g1(33,9)=(-1);
  g1(33,52)=1;
  g1(34,4)=(-1);
  g1(34,53)=1;
  g1(35,3)=(-1);
  g1(35,54)=1;
  g1(36,10)=(-1);
  g1(36,55)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],36,4761);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],36,328509);
end
end
end
end
