function [residual, g1, g2, g3] = SW_BLE_Estimation_Mode_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(38, 1);
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
forecast_i__ = params(35)^2*y(4)+(1+params(35))*params(36)*y(19);
forecast_q__ = params(37)^2*y(2)+(1+params(37))*params(38)*y(20);
forecast_rk__ = params(39)^2*y(1)+(1+params(39))*params(40)*y(21);
forecast_pinf__ = params(41)^2*y(7)+(1+params(41))*params(42)*y(22);
forecast_c__ = params(43)^2*y(3)+(1+params(43))*params(44)*y(23);
forecast_l__ = params(45)^2*y(6)+(1+params(45))*params(46)*y(24);
forecast_w__ = params(47)^2*y(8)+(1+params(47))*params(48)*y(25);
T171 = 1/(1+gamma__*beta_bar__);
T175 = gamma__^2;
T187 = (1-params(4))/(1-params(4)+Rk__);
T199 = params(13)/gamma__;
T203 = (1-T199)/(params(19)*(1+T199));
T219 = (params(19)-1)*WL_C__/(params(19)*(1+T199));
T250 = 1/(1+gamma__*beta_bar__*params(17));
T264 = (1-params(18))*(1-gamma__*beta_bar__*params(18))/params(18)/(1+(params(14)-1)*params(2));
T272 = gamma__*beta_bar__/(1+gamma__*beta_bar__);
T299 = (1-params(16))*(1-gamma__*beta_bar__*params(16))/((1+gamma__*beta_bar__)*params(16))*1/(1+(params(5)-1)*params(1));
lhs =y(35);
rhs =params(9)*y(37)+(1-params(9))*y(45)-y(48);
residual(1)= lhs-rhs;
lhs =y(36);
rhs =y(37)*(1-params(11))/params(11);
residual(2)= lhs-rhs;
lhs =y(37);
rhs =y(45)+y(43)-y(38);
residual(3)= lhs-rhs;
lhs =y(38);
rhs =y(36)+y(18);
residual(4)= lhs-rhs;
lhs =y(41);
rhs =T171*(y(4)+gamma__*beta_bar__*forecast_i__+1/(T175*params(12))*y(39))+y(51);
residual(5)= lhs-rhs;
lhs =y(39);
rhs =T187*forecast_q__+Rk__/(1-params(4)+Rk__)*forecast_rk__-y(46)+forecast_pinf__+1/T203*y(49);
residual(6)= lhs-rhs;
lhs =y(40);
rhs =y(49)+y(3)*T199/(1+T199)+1/(1+T199)*forecast_c__+T219*(y(43)-forecast_l__)-T203*(y(46)-forecast_pinf__);
residual(7)= lhs-rhs;
lhs =y(42);
rhs =C_Y__*y(40)+y(41)*I_Y__+y(50)+y(36)*Z_Y__;
residual(8)= lhs-rhs;
lhs =y(42);
rhs =params(14)*(y(48)+params(9)*y(38)+(1-params(9))*y(43));
residual(9)= lhs-rhs;
lhs =y(44);
rhs =T250*(gamma__*beta_bar__*forecast_pinf__+y(7)*params(17)+y(35)*T264)+y(53);
residual(10)= lhs-rhs;
lhs =y(45);
rhs =y(8)*T171+T272*forecast_w__+y(7)*params(15)/(1+gamma__*beta_bar__)-y(44)*(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__)+forecast_pinf__*T272+T299*(y(43)*params(20)+y(40)*1/(1-T199)-y(3)*T199/(1-T199)-y(45))+y(54);
residual(11)= lhs-rhs;
lhs =y(46);
rhs =y(44)*params(21)*(1-params(24))+(1-params(24))*params(23)*(y(42)-params(14)*y(48))+params(22)*(y(42)-params(14)*y(48)-(y(5)-params(14)*y(11)))+params(24)*y(9)+y(52);
residual(12)= lhs-rhs;
lhs =y(48);
rhs =y(11)*params(25)+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(50);
rhs =params(32)*y(13)+x(it_, 3)+x(it_, 1)*params(31);
residual(14)= lhs-rhs;
lhs =y(51);
rhs =params(29)*y(14)+x(it_, 4);
residual(15)= lhs-rhs;
lhs =y(49);
rhs =params(26)*y(12)+x(it_, 2);
residual(16)= lhs-rhs;
lhs =y(53);
rhs =params(27)*y(16)+x(it_, 6)-params(34)*y(26);
residual(17)= lhs-rhs;
lhs =y(54);
rhs =params(28)*y(17)+x(it_, 7)-params(33)*y(27);
residual(18)= lhs-rhs;
lhs =y(52);
rhs =params(30)*y(15)+x(it_, 5);
residual(19)= lhs-rhs;
lhs =y(55);
rhs =y(18)*(1-I_K_bar__)+y(41)*I_K_bar__+y(51)*params(12)*T175*I_K_bar__;
residual(20)= lhs-rhs;
lhs =y(47);
rhs =y(42)-params(14)*y(48);
residual(21)= lhs-rhs;
lhs =y(31);
rhs =params(10)+y(42)-y(5);
residual(22)= lhs-rhs;
lhs =y(32);
rhs =params(10)+y(40)-y(3);
residual(23)= lhs-rhs;
lhs =y(33);
rhs =params(10)+y(41)-y(4);
residual(24)= lhs-rhs;
lhs =y(34);
rhs =params(10)+y(45)-y(8);
residual(25)= lhs-rhs;
lhs =y(30);
rhs =params(7)+y(44);
residual(26)= lhs-rhs;
lhs =y(29);
rhs =y(46)+r_bar__;
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(43)+params(6);
residual(28)= lhs-rhs;
lhs =y(56);
rhs =y(47)-y(10);
residual(29)= lhs-rhs;
lhs =y(57);
rhs =y(4);
residual(30)= lhs-rhs;
lhs =y(58);
rhs =y(2);
residual(31)= lhs-rhs;
lhs =y(59);
rhs =y(1);
residual(32)= lhs-rhs;
lhs =y(60);
rhs =y(7);
residual(33)= lhs-rhs;
lhs =y(61);
rhs =y(3);
residual(34)= lhs-rhs;
lhs =y(62);
rhs =y(6);
residual(35)= lhs-rhs;
lhs =y(63);
rhs =y(8);
residual(36)= lhs-rhs;
lhs =y(64);
rhs =x(it_, 6);
residual(37)= lhs-rhs;
lhs =y(65);
rhs =x(it_, 7);
residual(38)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(38, 72);

  %
  % Jacobian matrix
  %

T105 = params(41)^2;
T521 = (-T105);
  g1(1,35)=1;
  g1(1,37)=(-params(9));
  g1(1,45)=(-(1-params(9)));
  g1(1,48)=1;
  g1(2,36)=1;
  g1(2,37)=(-((1-params(11))/params(11)));
  g1(3,37)=1;
  g1(3,38)=1;
  g1(3,43)=(-1);
  g1(3,45)=(-1);
  g1(4,36)=(-1);
  g1(4,38)=1;
  g1(4,18)=(-1);
  g1(5,39)=(-(T171*1/(T175*params(12))));
  g1(5,4)=(-(T171*(1+params(35)^2*gamma__*beta_bar__)));
  g1(5,41)=1;
  g1(5,51)=(-1);
  g1(5,19)=(-(T171*(1+params(35))*params(36)*gamma__*beta_bar__));
  g1(6,1)=(-(params(39)^2*Rk__/(1-params(4)+Rk__)));
  g1(6,2)=(-(params(37)^2*T187));
  g1(6,39)=1;
  g1(6,7)=T521;
  g1(6,46)=1;
  g1(6,49)=(-(1/T203));
  g1(6,20)=(-((1+params(37))*params(38)*T187));
  g1(6,21)=(-((1+params(39))*params(40)*Rk__/(1-params(4)+Rk__)));
  g1(6,22)=(-((1+params(41))*params(42)));
  g1(7,3)=(-(T199/(1+T199)+params(43)^2*1/(1+T199)));
  g1(7,40)=1;
  g1(7,6)=(-(T219*(-(params(45)^2))));
  g1(7,43)=(-T219);
  g1(7,7)=T203*T521;
  g1(7,46)=T203;
  g1(7,49)=(-1);
  g1(7,22)=T203*(-((1+params(41))*params(42)));
  g1(7,23)=(-((1+params(43))*params(44)*1/(1+T199)));
  g1(7,24)=(-(T219*(-((1+params(45))*params(46)))));
  g1(8,36)=(-Z_Y__);
  g1(8,40)=(-C_Y__);
  g1(8,41)=(-I_Y__);
  g1(8,42)=1;
  g1(8,50)=(-1);
  g1(9,38)=(-(params(9)*params(14)));
  g1(9,42)=1;
  g1(9,43)=(-((1-params(9))*params(14)));
  g1(9,48)=(-params(14));
  g1(10,35)=(-(T250*T264));
  g1(10,7)=(-(T250*(params(17)+T105*gamma__*beta_bar__)));
  g1(10,44)=1;
  g1(10,53)=(-1);
  g1(10,22)=(-(T250*(1+params(41))*params(42)*gamma__*beta_bar__));
  g1(11,3)=(-(T299*(-(T199/(1-T199)))));
  g1(11,40)=(-(T299*1/(1-T199)));
  g1(11,43)=(-(T299*params(20)));
  g1(11,7)=(-(params(15)/(1+gamma__*beta_bar__)+T105*T272));
  g1(11,44)=(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__);
  g1(11,8)=(-(T171+params(47)^2*T272));
  g1(11,45)=1-(-T299);
  g1(11,54)=(-1);
  g1(11,22)=(-((1+params(41))*params(42)*T272));
  g1(11,25)=(-((1+params(47))*params(48)*T272));
  g1(12,5)=params(22);
  g1(12,42)=(-((1-params(24))*params(23)+params(22)));
  g1(12,44)=(-(params(21)*(1-params(24))));
  g1(12,9)=(-params(24));
  g1(12,46)=1;
  g1(12,11)=(-(params(14)*params(22)));
  g1(12,48)=(-((1-params(24))*params(23)*(-params(14))+params(22)*(-params(14))));
  g1(12,52)=(-1);
  g1(13,11)=(-params(25));
  g1(13,48)=1;
  g1(13,66)=(-1);
  g1(14,13)=(-params(32));
  g1(14,50)=1;
  g1(14,66)=(-params(31));
  g1(14,68)=(-1);
  g1(15,14)=(-params(29));
  g1(15,51)=1;
  g1(15,69)=(-1);
  g1(16,12)=(-params(26));
  g1(16,49)=1;
  g1(16,67)=(-1);
  g1(17,16)=(-params(27));
  g1(17,53)=1;
  g1(17,71)=(-1);
  g1(17,26)=params(34);
  g1(18,17)=(-params(28));
  g1(18,54)=1;
  g1(18,72)=(-1);
  g1(18,27)=params(33);
  g1(19,15)=(-params(30));
  g1(19,52)=1;
  g1(19,70)=(-1);
  g1(20,41)=(-I_K_bar__);
  g1(20,51)=(-(params(12)*T175*I_K_bar__));
  g1(20,18)=(-(1-I_K_bar__));
  g1(20,55)=1;
  g1(21,42)=(-1);
  g1(21,47)=1;
  g1(21,48)=params(14);
  g1(22,31)=1;
  g1(22,5)=1;
  g1(22,42)=(-1);
  g1(23,32)=1;
  g1(23,3)=1;
  g1(23,40)=(-1);
  g1(24,33)=1;
  g1(24,4)=1;
  g1(24,41)=(-1);
  g1(25,34)=1;
  g1(25,8)=1;
  g1(25,45)=(-1);
  g1(26,30)=1;
  g1(26,44)=(-1);
  g1(27,29)=1;
  g1(27,46)=(-1);
  g1(28,28)=1;
  g1(28,43)=(-1);
  g1(29,10)=1;
  g1(29,47)=(-1);
  g1(29,56)=1;
  g1(30,4)=(-1);
  g1(30,57)=1;
  g1(31,2)=(-1);
  g1(31,58)=1;
  g1(32,1)=(-1);
  g1(32,59)=1;
  g1(33,7)=(-1);
  g1(33,60)=1;
  g1(34,3)=(-1);
  g1(34,61)=1;
  g1(35,6)=(-1);
  g1(35,62)=1;
  g1(36,8)=(-1);
  g1(36,63)=1;
  g1(37,71)=(-1);
  g1(37,64)=1;
  g1(38,72)=(-1);
  g1(38,65)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],38,5184);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],38,373248);
end
end
end
end
