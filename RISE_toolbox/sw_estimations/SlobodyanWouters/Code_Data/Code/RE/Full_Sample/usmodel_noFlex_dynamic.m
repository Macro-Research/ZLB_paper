function [residual, g1, g2, g3] = usmodel_noFlex_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(40, 1);
T21 = 1/(params(10)/(1-params(10)));
T37 = 1/(1+params(44)*params(42));
T42 = params(42)^2;
T44 = T42*params(12);
T57 = params(15)/params(42);
T62 = (1-T57)/(params(14)*(1+T57));
T74 = (1-params(13))/(params(47)+1-params(13));
T90 = (params(14)-1)*params(56)/(params(14)*(1+T57));
T120 = 1/(1-T57);
T123 = T57/(1-T57);
T212 = 1/(1+params(44)*params(42)*params(21));
T227 = (1-params(22))*(1-params(44)*params(42)*params(22))/params(22)/(1+(params(18)-1)*params(3));
T236 = params(44)*params(42)/(1+params(44)*params(42));
T264 = (1-params(20))*(1-params(44)*params(42)*params(20))/((1+params(44)*params(42))*params(20))*1/(1+(params(24)-1)*params(1));
lhs =y(51);
rhs =params(9)*y(30)+(1-params(9))*y(37);
residual(1)= lhs-rhs;
lhs =y(29);
rhs =y(30)*T21;
residual(2)= lhs-rhs;
lhs =y(30);
rhs =y(37)+y(36)-y(31);
residual(3)= lhs-rhs;
lhs =y(31);
rhs =y(29)+y(18);
residual(4)= lhs-rhs;
lhs =y(34);
rhs =T37*(y(4)+params(44)*params(42)*y(63)+1/T44*y(32))+y(54);
residual(5)= lhs-rhs;
lhs =y(32);
rhs =(-y(38))+y(52)*1/T62+params(47)/(params(47)+1-params(13))*y(60)+T74*y(61);
residual(6)= lhs-rhs;
lhs =y(33);
rhs =y(52)+T57/(1+T57)*y(3)+1/(1+T57)*y(62)+T90*(y(36)-y(64))-y(38)*T62;
residual(7)= lhs-rhs;
lhs =y(35);
rhs =y(33)*params(54)+y(34)*params(53)+y(53)+y(29)*params(55);
residual(8)= lhs-rhs;
lhs =y(35);
rhs =params(18)*(y(51)+params(9)*y(31)+(1-params(9))*y(36));
residual(9)= lhs-rhs;
lhs =y(37);
rhs =y(36)*params(23)+y(33)*T120-y(3)*T123;
residual(10)= lhs-rhs;
lhs =y(58);
rhs =y(18)*(1-params(49))+y(34)*params(49)+y(54)*T44*params(49);
residual(11)= lhs-rhs;
lhs =y(39);
rhs =params(9)*y(41)+(1-params(9))*y(49)-y(51);
residual(12)= lhs-rhs;
lhs =y(40);
rhs =T21*y(41);
residual(13)= lhs-rhs;
lhs =y(41);
rhs =y(49)+y(47)-y(42);
residual(14)= lhs-rhs;
lhs =y(42);
rhs =y(40)+y(19);
residual(15)= lhs-rhs;
lhs =y(45);
rhs =y(54)+T37*(y(6)+params(44)*params(42)*y(68)+1/T44*y(43));
residual(16)= lhs-rhs;
lhs =y(43);
rhs =y(52)*1/T62+(-y(50))+y(70)+params(47)/(params(47)+1-params(13))*y(65)+T74*y(66);
residual(17)= lhs-rhs;
lhs =y(44);
rhs =y(52)+T57/(1+T57)*y(5)+1/(1+T57)*y(67)+T90*(y(47)-y(69))-T62*(y(50)-y(70));
residual(18)= lhs-rhs;
lhs =y(46);
rhs =y(53)+params(54)*y(44)+params(53)*y(45)+params(55)*y(40);
residual(19)= lhs-rhs;
lhs =y(46);
rhs =params(18)*(y(51)+params(9)*y(42)+(1-params(9))*y(47));
residual(20)= lhs-rhs;
lhs =y(48);
rhs =T212*(params(44)*params(42)*y(70)+params(21)*y(8)+y(39)*T227)+y(56);
residual(21)= lhs-rhs;
lhs =y(49);
rhs =T37*y(9)+T236*y(71)+y(8)*params(19)/(1+params(44)*params(42))-y(48)*(1+params(44)*params(42)*params(19))/(1+params(44)*params(42))+y(70)*T236+T264*(params(23)*y(47)+T120*y(44)-T123*y(5)-y(49))+y(57);
residual(22)= lhs-rhs;
lhs =y(50);
rhs =y(48)*params(26)*(1-params(29))+(1-params(29))*params(28)*(y(46)-y(51)*params(18))+params(27)*(y(46)-y(51)*params(18)-y(7)+params(18)*y(11))+params(29)*y(10)+y(55);
residual(23)= lhs-rhs;
lhs =y(51);
rhs =y(11)*params(30)+x(it_, 1);
residual(24)= lhs-rhs;
lhs =y(52);
rhs =params(32)*y(12)+x(it_, 2);
residual(25)= lhs-rhs;
lhs =y(53);
rhs =params(33)*y(13)+x(it_, 3)+x(it_, 1)*params(2);
residual(26)= lhs-rhs;
lhs =y(54);
rhs =params(35)*y(14)+x(it_, 4);
residual(27)= lhs-rhs;
lhs =y(55);
rhs =params(36)*y(15)+x(it_, 5);
residual(28)= lhs-rhs;
lhs =y(56);
rhs =params(37)*y(16)+y(28)-params(8)*y(2);
residual(29)= lhs-rhs;
lhs =y(28);
rhs =x(it_, 6);
residual(30)= lhs-rhs;
lhs =y(57);
rhs =params(38)*y(17)+y(27)-params(7)*y(1);
residual(31)= lhs-rhs;
lhs =y(27);
rhs =x(it_, 7);
residual(32)= lhs-rhs;
lhs =y(59);
rhs =(1-params(49))*y(19)+params(49)*y(45)+y(54)*params(12)*T42*params(49);
residual(33)= lhs-rhs;
lhs =y(23);
rhs =y(46)-y(7)+params(39);
residual(34)= lhs-rhs;
lhs =y(24);
rhs =params(39)+y(44)-y(5);
residual(35)= lhs-rhs;
lhs =y(25);
rhs =params(39)+y(45)-y(6);
residual(36)= lhs-rhs;
lhs =y(26);
rhs =params(39)+y(49)-y(9);
residual(37)= lhs-rhs;
lhs =y(22);
rhs =y(48)+params(5);
residual(38)= lhs-rhs;
lhs =y(21);
rhs =y(50)+params(40);
residual(39)= lhs-rhs;
lhs =y(20);
rhs =y(47)+params(4);
residual(40)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(40, 78);

  %
  % Jacobian matrix
  %

  g1(1,30)=(-params(9));
  g1(1,37)=(-(1-params(9)));
  g1(1,51)=1;
  g1(2,29)=1;
  g1(2,30)=(-T21);
  g1(3,30)=1;
  g1(3,31)=1;
  g1(3,36)=(-1);
  g1(3,37)=(-1);
  g1(4,29)=(-1);
  g1(4,31)=1;
  g1(4,18)=(-1);
  g1(5,32)=(-(T37*1/T44));
  g1(5,4)=(-T37);
  g1(5,34)=1;
  g1(5,63)=(-(params(44)*params(42)*T37));
  g1(5,54)=(-1);
  g1(6,60)=(-(params(47)/(params(47)+1-params(13))));
  g1(6,32)=1;
  g1(6,61)=(-T74);
  g1(6,38)=1;
  g1(6,52)=(-(1/T62));
  g1(7,3)=(-(T57/(1+T57)));
  g1(7,33)=1;
  g1(7,62)=(-(1/(1+T57)));
  g1(7,36)=(-T90);
  g1(7,64)=T90;
  g1(7,38)=T62;
  g1(7,52)=(-1);
  g1(8,29)=(-params(55));
  g1(8,33)=(-params(54));
  g1(8,34)=(-params(53));
  g1(8,35)=1;
  g1(8,53)=(-1);
  g1(9,31)=(-(params(9)*params(18)));
  g1(9,35)=1;
  g1(9,36)=(-((1-params(9))*params(18)));
  g1(9,51)=(-params(18));
  g1(10,3)=T123;
  g1(10,33)=(-T120);
  g1(10,36)=(-params(23));
  g1(10,37)=1;
  g1(11,34)=(-params(49));
  g1(11,54)=(-(T44*params(49)));
  g1(11,18)=(-(1-params(49)));
  g1(11,58)=1;
  g1(12,39)=1;
  g1(12,41)=(-params(9));
  g1(12,49)=(-(1-params(9)));
  g1(12,51)=1;
  g1(13,40)=1;
  g1(13,41)=(-T21);
  g1(14,41)=1;
  g1(14,42)=1;
  g1(14,47)=(-1);
  g1(14,49)=(-1);
  g1(15,40)=(-1);
  g1(15,42)=1;
  g1(15,19)=(-1);
  g1(16,43)=(-(T37*1/T44));
  g1(16,6)=(-T37);
  g1(16,45)=1;
  g1(16,68)=(-(params(44)*params(42)*T37));
  g1(16,54)=(-1);
  g1(17,65)=(-(params(47)/(params(47)+1-params(13))));
  g1(17,43)=1;
  g1(17,66)=(-T74);
  g1(17,70)=(-1);
  g1(17,50)=1;
  g1(17,52)=(-(1/T62));
  g1(18,5)=(-(T57/(1+T57)));
  g1(18,44)=1;
  g1(18,67)=(-(1/(1+T57)));
  g1(18,47)=(-T90);
  g1(18,69)=T90;
  g1(18,70)=(-T62);
  g1(18,50)=T62;
  g1(18,52)=(-1);
  g1(19,40)=(-params(55));
  g1(19,44)=(-params(54));
  g1(19,45)=(-params(53));
  g1(19,46)=1;
  g1(19,53)=(-1);
  g1(20,42)=(-(params(9)*params(18)));
  g1(20,46)=1;
  g1(20,47)=(-((1-params(9))*params(18)));
  g1(20,51)=(-params(18));
  g1(21,39)=(-(T212*T227));
  g1(21,8)=(-(params(21)*T212));
  g1(21,48)=1;
  g1(21,70)=(-(params(44)*params(42)*T212));
  g1(21,56)=(-1);
  g1(22,5)=(-(T264*(-T123)));
  g1(22,44)=(-(T120*T264));
  g1(22,47)=(-(params(23)*T264));
  g1(22,8)=(-(params(19)/(1+params(44)*params(42))));
  g1(22,48)=(1+params(44)*params(42)*params(19))/(1+params(44)*params(42));
  g1(22,70)=(-T236);
  g1(22,9)=(-T37);
  g1(22,49)=1-(-T264);
  g1(22,71)=(-T236);
  g1(22,57)=(-1);
  g1(23,7)=params(27);
  g1(23,46)=(-((1-params(29))*params(28)+params(27)));
  g1(23,48)=(-(params(26)*(1-params(29))));
  g1(23,10)=(-params(29));
  g1(23,50)=1;
  g1(23,11)=(-(params(18)*params(27)));
  g1(23,51)=(-((1-params(29))*params(28)*(-params(18))+params(27)*(-params(18))));
  g1(23,55)=(-1);
  g1(24,11)=(-params(30));
  g1(24,51)=1;
  g1(24,72)=(-1);
  g1(25,12)=(-params(32));
  g1(25,52)=1;
  g1(25,73)=(-1);
  g1(26,13)=(-params(33));
  g1(26,53)=1;
  g1(26,72)=(-params(2));
  g1(26,74)=(-1);
  g1(27,14)=(-params(35));
  g1(27,54)=1;
  g1(27,75)=(-1);
  g1(28,15)=(-params(36));
  g1(28,55)=1;
  g1(28,76)=(-1);
  g1(29,2)=params(8);
  g1(29,28)=(-1);
  g1(29,16)=(-params(37));
  g1(29,56)=1;
  g1(30,28)=1;
  g1(30,77)=(-1);
  g1(31,1)=params(7);
  g1(31,27)=(-1);
  g1(31,17)=(-params(38));
  g1(31,57)=1;
  g1(32,27)=1;
  g1(32,78)=(-1);
  g1(33,45)=(-params(49));
  g1(33,54)=(-(params(12)*T42*params(49)));
  g1(33,19)=(-(1-params(49)));
  g1(33,59)=1;
  g1(34,23)=1;
  g1(34,7)=1;
  g1(34,46)=(-1);
  g1(35,24)=1;
  g1(35,5)=1;
  g1(35,44)=(-1);
  g1(36,25)=1;
  g1(36,6)=1;
  g1(36,45)=(-1);
  g1(37,26)=1;
  g1(37,9)=1;
  g1(37,49)=(-1);
  g1(38,22)=1;
  g1(38,48)=(-1);
  g1(39,21)=1;
  g1(39,50)=(-1);
  g1(40,20)=1;
  g1(40,47)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],40,6084);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],40,474552);
end
end
end
end
