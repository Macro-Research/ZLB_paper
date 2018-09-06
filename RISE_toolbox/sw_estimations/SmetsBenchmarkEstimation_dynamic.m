function [residual, g1, g2, g3] = SmetsBenchmarkEstimation_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(30, 1);
PI_star__ = 1+params(7)/100;
gamma__ = 1+params(10)/100;
beta__ = 1/(1+params(8)/100);
beta_bar__ = beta__*gamma__^(-params(19));
Rk__ = beta__^(-1)*gamma__^params(19)-(1-params(4));
W__ = (params(9)^params(9)*(1-params(9))^(1-params(9))/(params(14)*Rk__^params(9)))^(1/(1-params(9)));
I_K_bar__ = 1-(1-params(4))/gamma__;
I_K__ = gamma__*(1-(1-params(4))/gamma__);
L_K__ = (1-params(9))/params(9)*Rk__/W__;
K_Y__ = params(14)*L_K__^(params(9)-1);
I_Y__ = I_K__*K_Y__;
C_Y__ = 1-params(3)-I_K__*K_Y__;
Z_Y__ = Rk__*K_Y__;
WL_C__ = K_Y__*Rk__*(1-params(9))*1/params(5)/params(9)/C_Y__;
r_bar__ = 100*(PI_star__/(beta__*gamma__^(-params(19)))-1);
T100 = 1/(1+gamma__*beta_bar__);
T105 = gamma__^2;
T129 = params(13)/gamma__;
T133 = (1-T129)/(params(19)*(1+T129));
T150 = (params(19)-1)*WL_C__/(params(19)*(1+T129));
T181 = 1/(1+gamma__*beta_bar__*params(17));
T196 = (1-params(18))*(1-gamma__*beta_bar__*params(18))/params(18)/(1+(params(14)-1)*params(2));
T205 = gamma__*beta_bar__/(1+gamma__*beta_bar__);
T232 = (1-params(16))*(1-gamma__*beta_bar__*params(16))/((1+gamma__*beta_bar__)*params(16))*1/(1+(params(5)-1)*params(1));
lhs =y(24);
rhs =params(9)*y(26)+(1-params(9))*y(34)-y(37);
residual(1)= lhs-rhs;
lhs =y(25);
rhs =y(26)*(1-params(11))/params(11);
residual(2)= lhs-rhs;
lhs =y(26);
rhs =y(34)+y(32)-y(27);
residual(3)= lhs-rhs;
lhs =y(27);
rhs =y(25)+y(14);
residual(4)= lhs-rhs;
lhs =y(30);
rhs =T100*(y(2)+gamma__*beta_bar__*y(50)+1/(T105*params(12))*y(28))+y(40);
residual(5)= lhs-rhs;
lhs =y(28);
rhs =(1-params(4))/(1-params(4)+Rk__)*y(48)+Rk__/(1-params(4)+Rk__)*y(47)-y(35)+y(52)+1/T133*y(38);
residual(6)= lhs-rhs;
lhs =y(29);
rhs =y(38)+T129/(1+T129)*y(1)+1/(1+T129)*y(49)+T150*(y(32)-y(51))-T133*(y(35)-y(52));
residual(7)= lhs-rhs;
lhs =y(31);
rhs =C_Y__*y(29)+y(30)*I_Y__+y(39)+y(25)*Z_Y__;
residual(8)= lhs-rhs;
lhs =y(31);
rhs =params(14)*(y(37)+params(9)*y(27)+(1-params(9))*y(32));
residual(9)= lhs-rhs;
lhs =y(33);
rhs =T181*(gamma__*beta_bar__*y(52)+params(17)*y(4)+y(24)*T196)+y(42);
residual(10)= lhs-rhs;
lhs =y(34);
rhs =T100*y(5)+T205*y(53)+y(4)*params(15)/(1+gamma__*beta_bar__)-y(33)*(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__)+y(52)*T205+T232*(y(32)*params(20)+y(29)*1/(1-T129)-y(1)*T129/(1-T129)-y(34))+y(43);
residual(11)= lhs-rhs;
lhs =y(35);
rhs =y(33)*params(21)*(1-params(24))+(1-params(24))*params(23)*(y(31)-params(14)*y(37))+params(22)*(y(31)-params(14)*y(37)-(y(3)-params(14)*y(7)))+params(24)*y(6)+y(41);
residual(12)= lhs-rhs;
lhs =y(37);
rhs =y(7)*params(25)+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(39);
rhs =params(32)*y(9)+x(it_, 3)+x(it_, 1)*params(31);
residual(14)= lhs-rhs;
lhs =y(40);
rhs =params(29)*y(10)+x(it_, 4);
residual(15)= lhs-rhs;
lhs =y(38);
rhs =params(26)*y(8)+x(it_, 2);
residual(16)= lhs-rhs;
lhs =y(42);
rhs =params(27)*y(12)+x(it_, 6)-params(34)*y(15);
residual(17)= lhs-rhs;
lhs =y(43);
rhs =params(28)*y(13)+x(it_, 7)-params(33)*y(16);
residual(18)= lhs-rhs;
lhs =y(41);
rhs =params(30)*y(11)+x(it_, 5);
residual(19)= lhs-rhs;
lhs =y(44);
rhs =y(14)*(1-I_K_bar__)+y(30)*I_K_bar__+y(40)*params(12)*T105*I_K_bar__;
residual(20)= lhs-rhs;
lhs =y(36);
rhs =y(31)-params(14)*y(37);
residual(21)= lhs-rhs;
lhs =y(20);
rhs =params(10)+y(31)-y(3);
residual(22)= lhs-rhs;
lhs =y(21);
rhs =params(10)+y(29)-y(1);
residual(23)= lhs-rhs;
lhs =y(22);
rhs =params(10)+y(30)-y(2);
residual(24)= lhs-rhs;
lhs =y(23);
rhs =params(10)+y(34)-y(5);
residual(25)= lhs-rhs;
lhs =y(19);
rhs =params(7)+y(33);
residual(26)= lhs-rhs;
lhs =y(18);
rhs =y(35)+r_bar__;
residual(27)= lhs-rhs;
lhs =y(17);
rhs =y(32)+params(6);
residual(28)= lhs-rhs;
lhs =y(45);
rhs =x(it_, 6);
residual(29)= lhs-rhs;
lhs =y(46);
rhs =x(it_, 7);
residual(30)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(30, 60);

  %
  % Jacobian matrix
  %

  g1(1,24)=1;
  g1(1,26)=(-params(9));
  g1(1,34)=(-(1-params(9)));
  g1(1,37)=1;
  g1(2,25)=1;
  g1(2,26)=(-((1-params(11))/params(11)));
  g1(3,26)=1;
  g1(3,27)=1;
  g1(3,32)=(-1);
  g1(3,34)=(-1);
  g1(4,25)=(-1);
  g1(4,27)=1;
  g1(4,14)=(-1);
  g1(5,28)=(-(T100*1/(T105*params(12))));
  g1(5,2)=(-T100);
  g1(5,30)=1;
  g1(5,50)=(-(gamma__*beta_bar__*T100));
  g1(5,40)=(-1);
  g1(6,47)=(-(Rk__/(1-params(4)+Rk__)));
  g1(6,28)=1;
  g1(6,48)=(-((1-params(4))/(1-params(4)+Rk__)));
  g1(6,52)=(-1);
  g1(6,35)=1;
  g1(6,38)=(-(1/T133));
  g1(7,1)=(-(T129/(1+T129)));
  g1(7,29)=1;
  g1(7,49)=(-(1/(1+T129)));
  g1(7,32)=(-T150);
  g1(7,51)=T150;
  g1(7,52)=(-T133);
  g1(7,35)=T133;
  g1(7,38)=(-1);
  g1(8,25)=(-Z_Y__);
  g1(8,29)=(-C_Y__);
  g1(8,30)=(-I_Y__);
  g1(8,31)=1;
  g1(8,39)=(-1);
  g1(9,27)=(-(params(9)*params(14)));
  g1(9,31)=1;
  g1(9,32)=(-((1-params(9))*params(14)));
  g1(9,37)=(-params(14));
  g1(10,24)=(-(T181*T196));
  g1(10,4)=(-(params(17)*T181));
  g1(10,33)=1;
  g1(10,52)=(-(gamma__*beta_bar__*T181));
  g1(10,42)=(-1);
  g1(11,1)=(-(T232*(-(T129/(1-T129)))));
  g1(11,29)=(-(T232*1/(1-T129)));
  g1(11,32)=(-(T232*params(20)));
  g1(11,4)=(-(params(15)/(1+gamma__*beta_bar__)));
  g1(11,33)=(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__);
  g1(11,52)=(-T205);
  g1(11,5)=(-T100);
  g1(11,34)=1-(-T232);
  g1(11,53)=(-T205);
  g1(11,43)=(-1);
  g1(12,3)=params(22);
  g1(12,31)=(-((1-params(24))*params(23)+params(22)));
  g1(12,33)=(-(params(21)*(1-params(24))));
  g1(12,6)=(-params(24));
  g1(12,35)=1;
  g1(12,7)=(-(params(14)*params(22)));
  g1(12,37)=(-((1-params(24))*params(23)*(-params(14))+params(22)*(-params(14))));
  g1(12,41)=(-1);
  g1(13,7)=(-params(25));
  g1(13,37)=1;
  g1(13,54)=(-1);
  g1(14,9)=(-params(32));
  g1(14,39)=1;
  g1(14,54)=(-params(31));
  g1(14,56)=(-1);
  g1(15,10)=(-params(29));
  g1(15,40)=1;
  g1(15,57)=(-1);
  g1(16,8)=(-params(26));
  g1(16,38)=1;
  g1(16,55)=(-1);
  g1(17,12)=(-params(27));
  g1(17,42)=1;
  g1(17,59)=(-1);
  g1(17,15)=params(34);
  g1(18,13)=(-params(28));
  g1(18,43)=1;
  g1(18,60)=(-1);
  g1(18,16)=params(33);
  g1(19,11)=(-params(30));
  g1(19,41)=1;
  g1(19,58)=(-1);
  g1(20,30)=(-I_K_bar__);
  g1(20,40)=(-(params(12)*T105*I_K_bar__));
  g1(20,14)=(-(1-I_K_bar__));
  g1(20,44)=1;
  g1(21,31)=(-1);
  g1(21,36)=1;
  g1(21,37)=params(14);
  g1(22,20)=1;
  g1(22,3)=1;
  g1(22,31)=(-1);
  g1(23,21)=1;
  g1(23,1)=1;
  g1(23,29)=(-1);
  g1(24,22)=1;
  g1(24,2)=1;
  g1(24,30)=(-1);
  g1(25,23)=1;
  g1(25,5)=1;
  g1(25,34)=(-1);
  g1(26,19)=1;
  g1(26,33)=(-1);
  g1(27,18)=1;
  g1(27,35)=(-1);
  g1(28,17)=1;
  g1(28,32)=(-1);
  g1(29,59)=(-1);
  g1(29,45)=1;
  g1(30,60)=(-1);
  g1(30,46)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],30,3600);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],30,216000);
end
end
end
end
