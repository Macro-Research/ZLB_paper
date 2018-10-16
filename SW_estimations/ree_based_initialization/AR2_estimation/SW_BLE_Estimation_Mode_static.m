function [residual, g1, g2, g3] = SW_BLE_Estimation_Mode_static(y, x, params)
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

residual = zeros( 38, 1);

%
% Model equations
%

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
forecast_i__ = params(35)^2*y(14)+(1+params(35))*params(36)*y(30);
forecast_q__ = params(37)^2*y(12)+(1+params(37))*params(38)*y(31);
forecast_rk__ = params(39)^2*y(10)+(1+params(39))*params(40)*y(32);
forecast_pinf__ = params(41)^2*y(17)+(1+params(41))*params(42)*y(33);
forecast_c__ = params(43)^2*y(13)+(1+params(43))*params(44)*y(34);
forecast_l__ = params(45)^2*y(16)+(1+params(45))*params(46)*y(35);
forecast_w__ = params(47)^2*y(18)+(1+params(47))*params(48)*y(36);
T166 = 1/(1+gamma__*beta_bar__);
T170 = gamma__^2;
T181 = (1-params(4))/(1-params(4)+Rk__);
T193 = params(13)/gamma__;
T197 = (1-T193)/(params(19)*(1+T193));
T212 = (params(19)-1)*WL_C__/(params(19)*(1+T193));
T242 = 1/(1+gamma__*beta_bar__*params(17));
T256 = (1-params(18))*(1-gamma__*beta_bar__*params(18))/params(18)/(1+(params(14)-1)*params(2));
T264 = gamma__*beta_bar__/(1+gamma__*beta_bar__);
T291 = (1-params(16))*(1-gamma__*beta_bar__*params(16))/((1+gamma__*beta_bar__)*params(16))*1/(1+(params(5)-1)*params(1));
lhs =y(8);
rhs =params(9)*y(10)+(1-params(9))*y(18)-y(21);
residual(1)= lhs-rhs;
lhs =y(9);
rhs =y(10)*(1-params(11))/params(11);
residual(2)= lhs-rhs;
lhs =y(10);
rhs =y(16)+y(18)-y(11);
residual(3)= lhs-rhs;
lhs =y(11);
rhs =y(9)+y(28);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =T166*(y(14)+gamma__*beta_bar__*forecast_i__+y(12)*1/(T170*params(12)))+y(24);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =T181*forecast_q__+Rk__/(1-params(4)+Rk__)*forecast_rk__-y(19)+forecast_pinf__+1/T197*y(22);
residual(6)= lhs-rhs;
lhs =y(13);
rhs =y(22)+y(13)*T193/(1+T193)+1/(1+T193)*forecast_c__+T212*(y(16)-forecast_l__)-T197*(y(19)-forecast_pinf__);
residual(7)= lhs-rhs;
lhs =y(15);
rhs =C_Y__*y(13)+y(14)*I_Y__+y(23)+y(9)*Z_Y__;
residual(8)= lhs-rhs;
lhs =y(15);
rhs =params(14)*(y(21)+params(9)*y(11)+(1-params(9))*y(16));
residual(9)= lhs-rhs;
lhs =y(17);
rhs =T242*(gamma__*beta_bar__*forecast_pinf__+y(17)*params(17)+y(8)*T256)+y(26);
residual(10)= lhs-rhs;
lhs =y(18);
rhs =y(18)*T166+T264*forecast_w__+y(17)*params(15)/(1+gamma__*beta_bar__)-y(17)*(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__)+forecast_pinf__*T264+T291*(y(16)*params(20)+y(13)*1/(1-T193)-y(13)*T193/(1-T193)-y(18))+y(27);
residual(11)= lhs-rhs;
lhs =y(19);
rhs =y(17)*params(21)*(1-params(24))+(1-params(24))*params(23)*(y(15)-params(14)*y(21))+y(19)*params(24)+y(25);
residual(12)= lhs-rhs;
lhs =y(21);
rhs =y(21)*params(25)+x(1);
residual(13)= lhs-rhs;
lhs =y(23);
rhs =y(23)*params(32)+x(3)+x(1)*params(31);
residual(14)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(29)+x(4);
residual(15)= lhs-rhs;
lhs =y(22);
rhs =y(22)*params(26)+x(2);
residual(16)= lhs-rhs;
lhs =y(26);
rhs =y(26)*params(27)+x(6)-x(6)*params(34);
residual(17)= lhs-rhs;
lhs =y(27);
rhs =y(27)*params(28)+x(7)-x(7)*params(33);
residual(18)= lhs-rhs;
lhs =y(25);
rhs =y(25)*params(30)+x(5);
residual(19)= lhs-rhs;
lhs =y(28);
rhs =y(28)*(1-I_K_bar__)+y(14)*I_K_bar__+y(24)*params(12)*T170*I_K_bar__;
residual(20)= lhs-rhs;
lhs =y(20);
rhs =y(15)-params(14)*y(21);
residual(21)= lhs-rhs;
lhs =y(4);
rhs =params(10);
residual(22)= lhs-rhs;
lhs =y(5);
rhs =params(10);
residual(23)= lhs-rhs;
lhs =y(6);
rhs =params(10);
residual(24)= lhs-rhs;
lhs =y(7);
rhs =params(10);
residual(25)= lhs-rhs;
lhs =y(3);
rhs =params(7)+y(17);
residual(26)= lhs-rhs;
lhs =y(2);
rhs =y(19)+r_bar__;
residual(27)= lhs-rhs;
lhs =y(1);
rhs =y(16)+params(6);
residual(28)= lhs-rhs;
residual(29) = y(29);
lhs =y(30);
rhs =y(14);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =y(12);
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(10);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =y(17);
residual(33)= lhs-rhs;
lhs =y(34);
rhs =y(13);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =y(16);
residual(35)= lhs-rhs;
lhs =y(36);
rhs =y(18);
residual(36)= lhs-rhs;
lhs =y(37);
rhs =x(6);
residual(37)= lhs-rhs;
lhs =y(38);
rhs =x(7);
residual(38)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(38, 38);

  %
  % Jacobian matrix
  %

T104 = params(41)^2;
T457 = (-T104);
  g1(1,8)=1;
  g1(1,10)=(-params(9));
  g1(1,18)=(-(1-params(9)));
  g1(1,21)=1;
  g1(2,9)=1;
  g1(2,10)=(-((1-params(11))/params(11)));
  g1(3,10)=1;
  g1(3,11)=1;
  g1(3,16)=(-1);
  g1(3,18)=(-1);
  g1(4,9)=(-1);
  g1(4,11)=1;
  g1(4,28)=(-1);
  g1(5,12)=(-(T166*1/(T170*params(12))));
  g1(5,14)=1-T166*(1+params(35)^2*gamma__*beta_bar__);
  g1(5,24)=(-1);
  g1(5,30)=(-(T166*(1+params(35))*params(36)*gamma__*beta_bar__));
  g1(6,10)=(-(params(39)^2*Rk__/(1-params(4)+Rk__)));
  g1(6,12)=1-params(37)^2*T181;
  g1(6,17)=T457;
  g1(6,19)=1;
  g1(6,22)=(-(1/T197));
  g1(6,31)=(-((1+params(37))*params(38)*T181));
  g1(6,32)=(-((1+params(39))*params(40)*Rk__/(1-params(4)+Rk__)));
  g1(6,33)=(-((1+params(41))*params(42)));
  g1(7,13)=1-(T193/(1+T193)+params(43)^2*1/(1+T193));
  g1(7,16)=(-(T212*(1-params(45)^2)));
  g1(7,17)=T197*T457;
  g1(7,19)=T197;
  g1(7,22)=(-1);
  g1(7,33)=T197*(-((1+params(41))*params(42)));
  g1(7,34)=(-((1+params(43))*params(44)*1/(1+T193)));
  g1(7,35)=(-(T212*(-((1+params(45))*params(46)))));
  g1(8,9)=(-Z_Y__);
  g1(8,13)=(-C_Y__);
  g1(8,14)=(-I_Y__);
  g1(8,15)=1;
  g1(8,23)=(-1);
  g1(9,11)=(-(params(9)*params(14)));
  g1(9,15)=1;
  g1(9,16)=(-((1-params(9))*params(14)));
  g1(9,21)=(-params(14));
  g1(10,8)=(-(T242*T256));
  g1(10,17)=1-T242*(params(17)+T104*gamma__*beta_bar__);
  g1(10,26)=(-1);
  g1(10,33)=(-(T242*(1+params(41))*params(42)*gamma__*beta_bar__));
  g1(11,13)=(-(T291*(1/(1-T193)-T193/(1-T193))));
  g1(11,16)=(-(T291*params(20)));
  g1(11,17)=(-(params(15)/(1+gamma__*beta_bar__)-(1+gamma__*beta_bar__*params(15))/(1+gamma__*beta_bar__)+T104*T264));
  g1(11,18)=1-(T166+params(47)^2*T264-T291);
  g1(11,27)=(-1);
  g1(11,33)=(-((1+params(41))*params(42)*T264));
  g1(11,36)=(-((1+params(47))*params(48)*T264));
  g1(12,15)=(-((1-params(24))*params(23)));
  g1(12,17)=(-(params(21)*(1-params(24))));
  g1(12,19)=1-params(24);
  g1(12,21)=(-((1-params(24))*params(23)*(-params(14))));
  g1(12,25)=(-1);
  g1(13,21)=1-params(25);
  g1(14,23)=1-params(32);
  g1(15,24)=1-params(29);
  g1(16,22)=1-params(26);
  g1(17,26)=1-params(27);
  g1(18,27)=1-params(28);
  g1(19,25)=1-params(30);
  g1(20,14)=(-I_K_bar__);
  g1(20,24)=(-(params(12)*T170*I_K_bar__));
  g1(20,28)=1-(1-I_K_bar__);
  g1(21,15)=(-1);
  g1(21,20)=1;
  g1(21,21)=params(14);
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
  g1(29,29)=1;
  g1(30,14)=(-1);
  g1(30,30)=1;
  g1(31,12)=(-1);
  g1(31,31)=1;
  g1(32,10)=(-1);
  g1(32,32)=1;
  g1(33,17)=(-1);
  g1(33,33)=1;
  g1(34,13)=(-1);
  g1(34,34)=1;
  g1(35,16)=(-1);
  g1(35,35)=1;
  g1(36,18)=(-1);
  g1(36,36)=1;
  g1(37,37)=1;
  g1(38,38)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],38,1444);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],38,54872);
end
end
end
end
