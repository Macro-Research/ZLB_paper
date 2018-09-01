function [residual, g1, g2, g3] = fs2000_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(15, 1);
T76 = (-(params(6)/(1-params(6))));
T87 = exp((params(3)+x(it_, 1))*(-params(1)));
T90 = y(3)^params(1);
T92 = y(13)^(-params(1));
T101 = y(13)^(1-params(1));
T102 = T90*T87*(1-params(1))*y(6)*params(2)*T101;
T146 = y(3)^(params(1)-1);
T152 = T101*params(1)*exp((-params(1))*(params(3)+log(y(8))))*T146+(1-params(7))*exp((-(params(3)+log(y(8)))));
T153 = y(6)*params(2)*T152;
lhs =y(18);
rhs =exp(params(3)+x(it_, 1));
residual(1)= lhs-rhs;
lhs =log(y(5));
rhs =(1-params(5))*log(params(4))+params(5)*log(y(1))+x(it_, 2);
residual(2)= lhs-rhs;
residual(3) = (-y(6))/(y(5)*y(21)*y(20))+y(22);
lhs =y(9);
rhs =y(14)/y(13);
residual(4)= lhs-rhs;
residual(5) = y(14)/y(13)+T76*y(6)*y(7)/(1-y(13));
lhs =y(10);
rhs =y(6)*(1-params(1))*T87*T90*T92/y(9);
residual(6)= lhs-rhs;
residual(7) = 1/(y(6)*y(7))-T102/(y(20)*y(21)*y(5)*y(14));
lhs =y(11)+y(7);
rhs =T101*T87*T90+y(3)*(1-params(7))*exp((-(params(3)+x(it_, 1))));
residual(8)= lhs-rhs;
lhs =y(6)*y(7);
rhs =y(5);
residual(9)= lhs-rhs;
lhs =y(5)-1+y(12);
rhs =y(14);
residual(10)= lhs-rhs;
lhs =y(8);
rhs =exp(x(it_, 1));
residual(11)= lhs-rhs;
lhs =y(17);
rhs =T87*T90*T101;
residual(12)= lhs-rhs;
lhs =y(15);
rhs =y(18)*y(17)/y(4);
residual(13)= lhs-rhs;
lhs =y(16);
rhs =y(1)*y(6)/y(2)/y(18);
residual(14)= lhs-rhs;
lhs =y(19);
rhs =T153/(y(5)*y(21)*y(20));
residual(15)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(15, 24);

  %
  % Jacobian matrix
  %

T255 = getPowerDeriv(y(3),params(1),1);
T290 = getPowerDeriv(y(13),1-params(1),1);
  g1(1,18)=1;
  g1(1,23)=(-exp(params(3)+x(it_, 1)));
  g1(2,1)=(-(params(5)*1/y(1)));
  g1(2,5)=1/y(5);
  g1(2,24)=(-1);
  g1(3,5)=(-((-y(6))*y(21)*y(20)))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20));
  g1(3,6)=(-1)/(y(5)*y(21)*y(20));
  g1(3,20)=(-((-y(6))*y(5)*y(21)))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20));
  g1(3,21)=(-((-y(6))*y(5)*y(20)))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20));
  g1(3,22)=1;
  g1(4,9)=1;
  g1(4,13)=(-((-y(14))/(y(13)*y(13))));
  g1(4,14)=(-(1/y(13)));
  g1(5,6)=T76*y(7)/(1-y(13));
  g1(5,7)=T76*y(6)/(1-y(13));
  g1(5,13)=(-y(14))/(y(13)*y(13))+T76*y(6)*y(7)/((1-y(13))*(1-y(13)));
  g1(5,14)=1/y(13);
  g1(6,6)=(-(T92*T90*(1-params(1))*T87/y(9)));
  g1(6,9)=(-((-(y(6)*(1-params(1))*T87*T90*T92))/(y(9)*y(9))));
  g1(6,10)=1;
  g1(6,3)=(-(T92*y(6)*(1-params(1))*T87*T255/y(9)));
  g1(6,13)=(-(y(6)*(1-params(1))*T87*T90*getPowerDeriv(y(13),(-params(1)),1)/y(9)));
  g1(6,23)=(-(T92*T90*y(6)*(1-params(1))*(-params(1))*T87/y(9)));
  g1(7,5)=(-((-(T102*y(20)*y(21)*y(14)))/(y(20)*y(21)*y(5)*y(14)*y(20)*y(21)*y(5)*y(14))));
  g1(7,6)=(-y(7))/(y(6)*y(7)*y(6)*y(7))-T101*T90*T87*params(2)*(1-params(1))/(y(20)*y(21)*y(5)*y(14));
  g1(7,20)=(-((-(T102*y(21)*y(5)*y(14)))/(y(20)*y(21)*y(5)*y(14)*y(20)*y(21)*y(5)*y(14))));
  g1(7,7)=(-y(6))/(y(6)*y(7)*y(6)*y(7));
  g1(7,21)=(-((-(T102*y(20)*y(5)*y(14)))/(y(20)*y(21)*y(5)*y(14)*y(20)*y(21)*y(5)*y(14))));
  g1(7,3)=(-(T101*T87*(1-params(1))*y(6)*params(2)*T255/(y(20)*y(21)*y(5)*y(14))));
  g1(7,13)=(-(T90*T87*(1-params(1))*y(6)*params(2)*T290/(y(20)*y(21)*y(5)*y(14))));
  g1(7,14)=(-((-(T102*y(20)*y(5)*y(21)))/(y(20)*y(21)*y(5)*y(14)*y(20)*y(21)*y(5)*y(14))));
  g1(7,23)=(-(T101*T90*(1-params(1))*y(6)*params(2)*(-params(1))*T87/(y(20)*y(21)*y(5)*y(14))));
  g1(8,7)=1;
  g1(8,3)=(-((1-params(7))*exp((-(params(3)+x(it_, 1))))+T101*T87*T255));
  g1(8,11)=1;
  g1(8,13)=(-(T87*T90*T290));
  g1(8,23)=(-(T101*T90*(-params(1))*T87+y(3)*(1-params(7))*(-exp((-(params(3)+x(it_, 1)))))));
  g1(9,5)=(-1);
  g1(9,6)=y(7);
  g1(9,7)=y(6);
  g1(10,5)=1;
  g1(10,12)=1;
  g1(10,14)=(-1);
  g1(11,8)=1;
  g1(11,23)=(-exp(x(it_, 1)));
  g1(12,3)=(-(T87*T101*T255));
  g1(12,13)=(-(T87*T90*T290));
  g1(12,17)=1;
  g1(12,23)=(-(T90*T101*(-params(1))*T87));
  g1(13,15)=1;
  g1(13,4)=(-((-(y(18)*y(17)))/(y(4)*y(4))));
  g1(13,17)=(-(y(18)/y(4)));
  g1(13,18)=(-(y(17)/y(4)));
  g1(14,1)=(-(y(6)/y(2)/y(18)));
  g1(14,2)=(-(y(1)*(-y(6))/(y(2)*y(2))/y(18)));
  g1(14,6)=(-(y(1)*1/y(2)/y(18)));
  g1(14,16)=1;
  g1(14,18)=(-((-(y(1)*y(6)/y(2)))/(y(18)*y(18))));
  g1(15,5)=(-((-(y(21)*y(20)*T153))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20))));
  g1(15,6)=(-(params(2)*T152/(y(5)*y(21)*y(20))));
  g1(15,20)=(-((-(T153*y(5)*y(21)))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20))));
  g1(15,21)=(-((-(T153*y(5)*y(20)))/(y(5)*y(21)*y(20)*y(5)*y(21)*y(20))));
  g1(15,8)=(-(y(6)*params(2)*(T101*T146*params(1)*exp((-params(1))*(params(3)+log(y(8))))*(-params(1))*1/y(8)+(1-params(7))*exp((-(params(3)+log(y(8)))))*(-(1/y(8))))/(y(5)*y(21)*y(20))));
  g1(15,3)=(-(y(6)*params(2)*T101*params(1)*exp((-params(1))*(params(3)+log(y(8))))*getPowerDeriv(y(3),params(1)-1,1)/(y(5)*y(21)*y(20))));
  g1(15,13)=(-(y(6)*params(2)*params(1)*exp((-params(1))*(params(3)+log(y(8))))*T146*T290/(y(5)*y(21)*y(20))));
  g1(15,19)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],15,576);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],15,13824);
end
end
end
end
