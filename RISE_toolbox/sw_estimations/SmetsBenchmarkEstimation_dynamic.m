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

residual = zeros(28, 1);
T38 = 1/(1+params(4)*params(2));
T43 = params(2)^2;
T70 = params(28)/params(2);
T75 = (1-T70)/(params(34)*(1+T70));
T92 = (params(34)-1)*params(14)/(params(34)*(1+T70));
T125 = 1/(1+params(4)*params(2)*params(32));
T140 = (1-params(33))*(1-params(4)*params(2)*params(33))/params(33)/(1+(params(29)-1)*params(17));
T149 = params(4)*params(2)/(1+params(4)*params(2));
T177 = (1-params(31))*(1-params(4)*params(2)*params(31))/((1+params(4)*params(2))*params(31))*1/(1+(params(20)-1)*params(16));
lhs =y(22);
rhs =params(24)*y(24)+(1-params(24))*y(32)-y(35);
residual(1)= lhs-rhs;
lhs =y(23);
rhs =y(24)*(1-params(26))/params(26);
residual(2)= lhs-rhs;
lhs =y(24);
rhs =y(32)+y(30)-y(25);
residual(3)= lhs-rhs;
lhs =y(25);
rhs =y(23)+y(14);
residual(4)= lhs-rhs;
lhs =y(28);
rhs =T38*(y(2)+params(4)*params(2)*y(46)+1/(T43*params(27))*y(26))+y(38);
residual(5)= lhs-rhs;
lhs =y(26);
rhs =(1-params(19))/(1-params(19)+params(5))*y(44)+params(5)/(1-params(19)+params(5))*y(43)-y(33)+y(48)+1/T75*y(36);
residual(6)= lhs-rhs;
lhs =y(27);
rhs =y(36)+T70/(1+T70)*y(1)+1/(1+T70)*y(45)+T92*(y(30)-y(47))-T75*(y(33)-y(48));
residual(7)= lhs-rhs;
lhs =y(29);
rhs =y(27)*params(12)+y(28)*params(11)+y(37)+y(23)*params(13);
residual(8)= lhs-rhs;
lhs =y(29);
rhs =params(29)*(y(35)+params(24)*y(25)+(1-params(24))*y(30));
residual(9)= lhs-rhs;
lhs =y(31);
rhs =T125*(params(4)*params(2)*y(48)+params(32)*y(4)+y(22)*T140)+y(40);
residual(10)= lhs-rhs;
lhs =y(32);
rhs =T38*y(5)+T149*y(49)+y(4)*params(30)/(1+params(4)*params(2))-y(31)*(1+params(4)*params(2)*params(30))/(1+params(4)*params(2))+y(48)*T149+T177*(y(30)*params(35)+y(27)*1/(1-T70)-y(1)*T70/(1-T70)-y(32))+y(41);
residual(11)= lhs-rhs;
lhs =y(33);
rhs =y(31)*params(36)*(1-params(39))+(1-params(39))*params(38)*(y(29)-y(35)*params(29))+params(37)*(y(29)-y(35)*params(29)-(y(3)-params(29)*y(7)))+params(39)*y(6)+y(39);
residual(12)= lhs-rhs;
lhs =y(35);
rhs =y(7)*params(40)+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(37);
rhs =params(47)*y(9)+x(it_, 3)+x(it_, 1)*params(46);
residual(14)= lhs-rhs;
lhs =y(38);
rhs =params(44)*y(10)+x(it_, 4);
residual(15)= lhs-rhs;
lhs =y(36);
rhs =params(41)*y(8)+x(it_, 2);
residual(16)= lhs-rhs;
lhs =y(40);
rhs =params(42)*y(12)+x(it_, 6)-params(49)*x(it_-1, 6);
residual(17)= lhs-rhs;
lhs =y(41);
rhs =params(43)*y(13)+x(it_, 7)-params(48)*x(it_-1, 7);
residual(18)= lhs-rhs;
lhs =y(39);
rhs =params(45)*y(11)+x(it_, 5);
residual(19)= lhs-rhs;
lhs =y(42);
rhs =y(14)*(1-params(7))+y(28)*params(7)+y(38)*params(27)*T43*params(7);
residual(20)= lhs-rhs;
lhs =y(34);
rhs =y(29)-y(35)*params(29);
residual(21)= lhs-rhs;
lhs =y(18);
rhs =y(29)-y(3)+params(25);
residual(22)= lhs-rhs;
lhs =y(19);
rhs =params(25)+y(27)-y(1);
residual(23)= lhs-rhs;
lhs =y(20);
rhs =params(25)+y(28)-y(2);
residual(24)= lhs-rhs;
lhs =y(21);
rhs =params(25)+y(32)-y(5);
residual(25)= lhs-rhs;
lhs =y(17);
rhs =y(31)+params(22);
residual(26)= lhs-rhs;
lhs =y(16);
rhs =y(33)+params(15);
residual(27)= lhs-rhs;
lhs =y(15);
rhs =y(30)+params(21);
residual(28)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(28, 56);

  %
  % Jacobian matrix
  %

  g1(1,22)=1;
  g1(1,24)=(-params(24));
  g1(1,32)=(-(1-params(24)));
  g1(1,35)=1;
  g1(2,23)=1;
  g1(2,24)=(-((1-params(26))/params(26)));
  g1(3,24)=1;
  g1(3,25)=1;
  g1(3,30)=(-1);
  g1(3,32)=(-1);
  g1(4,23)=(-1);
  g1(4,25)=1;
  g1(4,14)=(-1);
  g1(5,26)=(-(T38*1/(T43*params(27))));
  g1(5,2)=(-T38);
  g1(5,28)=1;
  g1(5,46)=(-(params(4)*params(2)*T38));
  g1(5,38)=(-1);
  g1(6,43)=(-(params(5)/(1-params(19)+params(5))));
  g1(6,26)=1;
  g1(6,44)=(-((1-params(19))/(1-params(19)+params(5))));
  g1(6,48)=(-1);
  g1(6,33)=1;
  g1(6,36)=(-(1/T75));
  g1(7,1)=(-(T70/(1+T70)));
  g1(7,27)=1;
  g1(7,45)=(-(1/(1+T70)));
  g1(7,30)=(-T92);
  g1(7,47)=T92;
  g1(7,48)=(-T75);
  g1(7,33)=T75;
  g1(7,36)=(-1);
  g1(8,23)=(-params(13));
  g1(8,27)=(-params(12));
  g1(8,28)=(-params(11));
  g1(8,29)=1;
  g1(8,37)=(-1);
  g1(9,25)=(-(params(24)*params(29)));
  g1(9,29)=1;
  g1(9,30)=(-((1-params(24))*params(29)));
  g1(9,35)=(-params(29));
  g1(10,22)=(-(T125*T140));
  g1(10,4)=(-(params(32)*T125));
  g1(10,31)=1;
  g1(10,48)=(-(params(4)*params(2)*T125));
  g1(10,40)=(-1);
  g1(11,1)=(-(T177*(-(T70/(1-T70)))));
  g1(11,27)=(-(T177*1/(1-T70)));
  g1(11,30)=(-(T177*params(35)));
  g1(11,4)=(-(params(30)/(1+params(4)*params(2))));
  g1(11,31)=(1+params(4)*params(2)*params(30))/(1+params(4)*params(2));
  g1(11,48)=(-T149);
  g1(11,5)=(-T38);
  g1(11,32)=1-(-T177);
  g1(11,49)=(-T149);
  g1(11,41)=(-1);
  g1(12,3)=params(37);
  g1(12,29)=(-((1-params(39))*params(38)+params(37)));
  g1(12,31)=(-(params(36)*(1-params(39))));
  g1(12,6)=(-params(39));
  g1(12,33)=1;
  g1(12,7)=(-(params(29)*params(37)));
  g1(12,35)=(-((1-params(39))*params(38)*(-params(29))+params(37)*(-params(29))));
  g1(12,39)=(-1);
  g1(13,7)=(-params(40));
  g1(13,35)=1;
  g1(13,50)=(-1);
  g1(14,9)=(-params(47));
  g1(14,37)=1;
  g1(14,50)=(-params(46));
  g1(14,52)=(-1);
  g1(15,10)=(-params(44));
  g1(15,38)=1;
  g1(15,53)=(-1);
  g1(16,8)=(-params(41));
  g1(16,36)=1;
  g1(16,51)=(-1);
  g1(17,12)=(-params(42));
  g1(17,40)=1;
  g1(17,55)=params(49);
  g1(17,55)=(-1);
  g1(18,13)=(-params(43));
  g1(18,41)=1;
  g1(18,56)=params(48);
  g1(18,56)=(-1);
  g1(19,11)=(-params(45));
  g1(19,39)=1;
  g1(19,54)=(-1);
  g1(20,28)=(-params(7));
  g1(20,38)=(-(params(27)*T43*params(7)));
  g1(20,14)=(-(1-params(7)));
  g1(20,42)=1;
  g1(21,29)=(-1);
  g1(21,34)=1;
  g1(21,35)=params(29);
  g1(22,18)=1;
  g1(22,3)=1;
  g1(22,29)=(-1);
  g1(23,19)=1;
  g1(23,1)=1;
  g1(23,27)=(-1);
  g1(24,20)=1;
  g1(24,2)=1;
  g1(24,28)=(-1);
  g1(25,21)=1;
  g1(25,5)=1;
  g1(25,32)=(-1);
  g1(26,17)=1;
  g1(26,31)=(-1);
  g1(27,16)=1;
  g1(27,33)=(-1);
  g1(28,15)=1;
  g1(28,30)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],28,3136);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],28,175616);
end
end
end
end
