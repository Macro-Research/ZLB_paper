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

residual = zeros(31, 1);
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
forecast_i__ = params(35)^2*y(4)+params(37)*y(9)+params(38)*y(3)+y(4)*params(39)+params(40)*y(5)+params(41)*y(7)+params(42)*y(8);
forecast_q__ = params(43)^2*y(2)+y(9)*params(45)+y(3)*params(46)+y(4)*params(47)+y(5)*params(48)+y(7)*params(49)+y(8)*params(50);
forecast_rk__ = params(51)^2*y(1)+y(9)*params(53)+y(3)*params(54)+y(4)*params(55)+y(5)*params(56)+y(7)*params(57)+y(8)*params(58);
forecast_pinf__ = y(7)*params(59)^2+y(9)*params(61)+y(3)*params(62)+y(4)*params(63)+y(5)*params(64)+y(7)*params(65)+y(8)*params(66);
forecast_c__ = y(3)*params(67)^2+y(9)*params(69)+y(3)*params(70)+y(4)*params(71)+y(5)*params(72)+y(7)*params(73)+y(8)*params(74);
forecast_l__ = params(75)^2*y(6)+y(9)*params(77)+y(3)*params(78)+y(4)*params(79)+y(5)*params(80)+y(7)*params(81)+y(8)*params(82);
forecast_w__ = y(8)*params(83)^2+y(9)*params(85)+y(3)*params(86)+y(4)*params(87)+y(5)*params(88)+y(7)*params(89)+y(8)*params(90);
T256 = 1/(1+gamma__*beta_bar__);
T260 = gamma__^2;
T272 = (1-params(4))/(1-params(4)+Rk__);
T275 = Rk__/(1-params(4)+Rk__);
T284 = params(13)/gamma__;
T288 = (1-T284)/(params(19)*(1+T284));
T297 = 1/(1+T284);
T304 = (params(19)-1)*WL_C__/(params(19)*(1+T284));
T335 = 1/(1+gamma__*beta_bar__*params(17));
T349 = (1-params(18))*(1-gamma__*beta_bar__*params(18))/params(18)/(1+(params(14)-1)*params(2));
T357 = gamma__*beta_bar__/(1+gamma__*beta_bar__);
T384 = (1-params(16))*(1-gamma__*beta_bar__*params(16))/((1+gamma__*beta_bar__)*params(16))*1/(1+(params(5)-1)*params(1));
lhs =y(28);
rhs =params(9)*y(30)+(1-params(9))*y(38)-y(41);
residual(1)= lhs-rhs;
lhs =y(29);
rhs =y(30)*(1-params(11))/params(11);
residual(2)= lhs-rhs;
lhs =y(30);
rhs =y(38)+y(36)-y(31);
residual(3)= lhs-rhs;
lhs =y(31);
rhs =y(29)+y(18);
residual(4)= lhs-rhs;
lhs =y(34);
rhs =T256*(y(4)+gamma__*beta_bar__*forecast_i__+1/(T260*params(12))*y(32))+y(44);
residual(5)= lhs-rhs;
lhs =y(32);
rhs =T272*forecast_q__+T275*forecast_rk__-y(39)+forecast_pinf__+1/T288*y(42);
residual(6)= lhs-rhs;
lhs =y(33);
rhs =y(42)+y(3)*T284/(1+T284)+T297*forecast_c__+T304*(y(36)-forecast_l__)-T288*(y(39)-forecast_pinf__);
residual(7)= lhs-rhs;
lhs =y(35);
rhs =C_Y__*y(33)+y(34)*I_Y__+y(43)+y(29)*Z_Y__;
residual(8)= lhs-rhs;
lhs =y(35);
rhs =params(14)*(y(41)+params(9)*y(31)+(1-params(9))*y(36));
residual(9)= lhs-rhs;
lhs =y(37);
rhs =T335*(gamma__*beta_bar__*forecast_pinf__+y(7)*params(17)+y(28)*T349)+y(46);
residual(10)= lhs-rhs;
lhs =y(38);
rhs =y(8)*T256+T357*forecast_w__+y(7)*params(15)/(1+gamma__*beta_bar__)-y(37)*(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__)+forecast_pinf__*T357+T384*(y(36)*params(20)+y(33)*1/(1-T284)-y(3)*T284/(1-T284)-y(38))+y(47);
residual(11)= lhs-rhs;
lhs =y(39);
rhs =y(37)*params(21)*(1-params(24))+(1-params(24))*params(23)*(y(35)-params(14)*y(41))+params(22)*(y(35)-params(14)*y(41)-(y(5)-params(14)*y(11)))+y(9)*params(24)+y(45);
residual(12)= lhs-rhs;
lhs =y(41);
rhs =y(11)*params(25)+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(43);
rhs =params(32)*y(13)+x(it_, 3)+x(it_, 1)*params(31);
residual(14)= lhs-rhs;
lhs =y(44);
rhs =params(29)*y(14)+x(it_, 4);
residual(15)= lhs-rhs;
lhs =y(42);
rhs =params(26)*y(12)+x(it_, 2);
residual(16)= lhs-rhs;
lhs =y(46);
rhs =params(27)*y(16)+x(it_, 6)-params(34)*y(19);
residual(17)= lhs-rhs;
lhs =y(47);
rhs =params(28)*y(17)+x(it_, 7)-params(33)*y(20);
residual(18)= lhs-rhs;
lhs =y(45);
rhs =params(30)*y(15)+x(it_, 5);
residual(19)= lhs-rhs;
lhs =y(48);
rhs =y(18)*(1-I_K_bar__)+y(34)*I_K_bar__+y(44)*params(12)*T260*I_K_bar__;
residual(20)= lhs-rhs;
lhs =y(40);
rhs =y(35)-params(14)*y(41);
residual(21)= lhs-rhs;
lhs =y(24);
rhs =params(10)+y(35)-y(5);
residual(22)= lhs-rhs;
lhs =y(25);
rhs =params(10)+y(33)-y(3);
residual(23)= lhs-rhs;
lhs =y(26);
rhs =params(10)+y(34)-y(4);
residual(24)= lhs-rhs;
lhs =y(27);
rhs =params(10)+y(38)-y(8);
residual(25)= lhs-rhs;
lhs =y(23);
rhs =params(7)+y(37);
residual(26)= lhs-rhs;
lhs =y(22);
rhs =y(39)+r_bar__;
residual(27)= lhs-rhs;
lhs =y(21);
rhs =y(36)+params(6);
residual(28)= lhs-rhs;
lhs =y(49);
rhs =y(40)-y(10);
residual(29)= lhs-rhs;
lhs =y(50);
rhs =x(it_, 6);
residual(30)= lhs-rhs;
lhs =y(51);
rhs =x(it_, 7);
residual(31)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(31, 58);

  %
  % Jacobian matrix
  %

T641 = params(59)^2+params(65);
  g1(1,28)=1;
  g1(1,30)=(-params(9));
  g1(1,38)=(-(1-params(9)));
  g1(1,41)=1;
  g1(2,29)=1;
  g1(2,30)=(-((1-params(11))/params(11)));
  g1(3,30)=1;
  g1(3,31)=1;
  g1(3,36)=(-1);
  g1(3,38)=(-1);
  g1(4,29)=(-1);
  g1(4,31)=1;
  g1(4,18)=(-1);
  g1(5,32)=(-(T256*1/(T260*params(12))));
  g1(5,3)=(-(T256*params(38)*gamma__*beta_bar__));
  g1(5,4)=(-(T256*(1+gamma__*beta_bar__*(params(35)^2+params(39)))));
  g1(5,34)=1;
  g1(5,5)=(-(T256*params(40)*gamma__*beta_bar__));
  g1(5,7)=(-(T256*params(41)*gamma__*beta_bar__));
  g1(5,8)=(-(T256*params(42)*gamma__*beta_bar__));
  g1(5,9)=(-(T256*params(37)*gamma__*beta_bar__));
  g1(5,44)=(-1);
  g1(6,1)=(-(params(51)^2*T275));
  g1(6,2)=(-(params(43)^2*T272));
  g1(6,32)=1;
  g1(6,3)=(-(params(62)+params(46)*T272+params(54)*T275));
  g1(6,4)=(-(params(63)+params(47)*T272+params(55)*T275));
  g1(6,5)=(-(params(64)+params(48)*T272+params(56)*T275));
  g1(6,7)=(-(params(49)*T272+params(57)*T275+T641));
  g1(6,8)=(-(params(66)+params(50)*T272+params(58)*T275));
  g1(6,9)=(-(params(61)+params(45)*T272+params(53)*T275));
  g1(6,39)=1;
  g1(6,42)=(-(1/T288));
  g1(7,3)=(-(T284/(1+T284)+T297*(params(67)^2+params(70))+T304*(-params(78))-T288*(-params(62))));
  g1(7,33)=1;
  g1(7,4)=(-(params(71)*T297+T304*(-params(79))-T288*(-params(63))));
  g1(7,5)=(-(params(72)*T297+T304*(-params(80))-T288*(-params(64))));
  g1(7,6)=(-(T304*(-(params(75)^2))));
  g1(7,36)=(-T304);
  g1(7,7)=(-(params(73)*T297+T304*(-params(81))-T288*(-T641)));
  g1(7,8)=(-(params(74)*T297+T304*(-params(82))-T288*(-params(66))));
  g1(7,9)=(-(params(69)*T297+T304*(-params(77))-T288*(-params(61))));
  g1(7,39)=T288;
  g1(7,42)=(-1);
  g1(8,29)=(-Z_Y__);
  g1(8,33)=(-C_Y__);
  g1(8,34)=(-I_Y__);
  g1(8,35)=1;
  g1(8,43)=(-1);
  g1(9,31)=(-(params(9)*params(14)));
  g1(9,35)=1;
  g1(9,36)=(-((1-params(9))*params(14)));
  g1(9,41)=(-params(14));
  g1(10,28)=(-(T335*T349));
  g1(10,3)=(-(T335*params(62)*gamma__*beta_bar__));
  g1(10,4)=(-(T335*params(63)*gamma__*beta_bar__));
  g1(10,5)=(-(T335*params(64)*gamma__*beta_bar__));
  g1(10,7)=(-(T335*(params(17)+gamma__*beta_bar__*T641)));
  g1(10,37)=1;
  g1(10,8)=(-(T335*params(66)*gamma__*beta_bar__));
  g1(10,9)=(-(T335*params(61)*gamma__*beta_bar__));
  g1(10,46)=(-1);
  g1(11,3)=(-(params(86)*T357+params(62)*T357+T384*(-(T284/(1-T284)))));
  g1(11,33)=(-(T384*1/(1-T284)));
  g1(11,4)=(-(params(87)*T357+params(63)*T357));
  g1(11,5)=(-(params(88)*T357+params(64)*T357));
  g1(11,36)=(-(T384*params(20)));
  g1(11,7)=(-(params(15)/(1+gamma__*beta_bar__)+params(89)*T357+T357*T641));
  g1(11,37)=(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__);
  g1(11,8)=(-(T256+T357*(params(83)^2+params(90))+params(66)*T357));
  g1(11,38)=1-(-T384);
  g1(11,9)=(-(params(85)*T357+params(61)*T357));
  g1(11,47)=(-1);
  g1(12,5)=params(22);
  g1(12,35)=(-((1-params(24))*params(23)+params(22)));
  g1(12,37)=(-(params(21)*(1-params(24))));
  g1(12,9)=(-params(24));
  g1(12,39)=1;
  g1(12,11)=(-(params(14)*params(22)));
  g1(12,41)=(-((1-params(24))*params(23)*(-params(14))+params(22)*(-params(14))));
  g1(12,45)=(-1);
  g1(13,11)=(-params(25));
  g1(13,41)=1;
  g1(13,52)=(-1);
  g1(14,13)=(-params(32));
  g1(14,43)=1;
  g1(14,52)=(-params(31));
  g1(14,54)=(-1);
  g1(15,14)=(-params(29));
  g1(15,44)=1;
  g1(15,55)=(-1);
  g1(16,12)=(-params(26));
  g1(16,42)=1;
  g1(16,53)=(-1);
  g1(17,16)=(-params(27));
  g1(17,46)=1;
  g1(17,57)=(-1);
  g1(17,19)=params(34);
  g1(18,17)=(-params(28));
  g1(18,47)=1;
  g1(18,58)=(-1);
  g1(18,20)=params(33);
  g1(19,15)=(-params(30));
  g1(19,45)=1;
  g1(19,56)=(-1);
  g1(20,34)=(-I_K_bar__);
  g1(20,44)=(-(params(12)*T260*I_K_bar__));
  g1(20,18)=(-(1-I_K_bar__));
  g1(20,48)=1;
  g1(21,35)=(-1);
  g1(21,40)=1;
  g1(21,41)=params(14);
  g1(22,24)=1;
  g1(22,5)=1;
  g1(22,35)=(-1);
  g1(23,25)=1;
  g1(23,3)=1;
  g1(23,33)=(-1);
  g1(24,26)=1;
  g1(24,4)=1;
  g1(24,34)=(-1);
  g1(25,27)=1;
  g1(25,8)=1;
  g1(25,38)=(-1);
  g1(26,23)=1;
  g1(26,37)=(-1);
  g1(27,22)=1;
  g1(27,39)=(-1);
  g1(28,21)=1;
  g1(28,36)=(-1);
  g1(29,10)=1;
  g1(29,40)=(-1);
  g1(29,49)=1;
  g1(30,57)=(-1);
  g1(30,50)=1;
  g1(31,58)=(-1);
  g1(31,51)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],31,3364);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],31,195112);
end
end
end
end
