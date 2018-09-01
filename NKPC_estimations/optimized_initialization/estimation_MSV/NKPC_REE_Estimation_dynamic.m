function [residual, g1, g2, g3] = NKPC_REE_Estimation_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(8, 1);
beta__ = 0.99;
y_forecast__ = params(11)+params(14)*(params(13)+params(20)*y(1)+params(21)*y(7)+params(22)*y(8))+y(7)*params(15)*params(8)+y(8)*params(16)*params(9);
pi_forecast__ = params(12)+(params(13)+params(20)*y(1)+params(21)*y(7)+params(22)*y(8))*params(17)+y(7)*params(8)*params(18)+y(8)*params(9)*params(19);
T51 = 1/params(5);
lhs =y(4);
rhs =y(7)+y_forecast__-T51*(y(6)-pi_forecast__);
residual(1)= lhs-rhs;
lhs =y(5);
rhs =y(8)+pi_forecast__*beta__+y(4)*params(4);
residual(2)= lhs-rhs;
lhs =y(6);
rhs =y(1)*params(10)+(1-params(10))*(y(5)*params(6)+y(4)*params(7))+x(it_, 3);
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(8)*y(2)+x(it_, 1);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(9)*y(3)+x(it_, 2);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(4)+params(1);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(5)+params(2);
residual(7)= lhs-rhs;
lhs =y(11);
rhs =y(6)+params(3);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 14);

  %
  % Jacobian matrix
  %

  g1(1,4)=1;
  g1(1,1)=(-(params(14)*params(20)-T51*(-(params(20)*params(17)))));
  g1(1,6)=T51;
  g1(1,7)=(-(1+params(15)*params(8)+params(14)*params(21)-T51*(-(params(8)*params(18)+params(21)*params(17)))));
  g1(1,8)=(-(params(16)*params(9)+params(14)*params(22)-T51*(-(params(9)*params(19)+params(22)*params(17)))));
  g1(2,4)=(-params(4));
  g1(2,5)=1;
  g1(2,1)=(-(beta__*params(20)*params(17)));
  g1(2,7)=(-(beta__*(params(8)*params(18)+params(21)*params(17))));
  g1(2,8)=(-(1+beta__*(params(9)*params(19)+params(22)*params(17))));
  g1(3,4)=(-((1-params(10))*params(7)));
  g1(3,5)=(-((1-params(10))*params(6)));
  g1(3,1)=(-params(10));
  g1(3,6)=1;
  g1(3,14)=(-1);
  g1(4,2)=(-params(8));
  g1(4,7)=1;
  g1(4,12)=(-1);
  g1(5,3)=(-params(9));
  g1(5,8)=1;
  g1(5,13)=(-1);
  g1(6,4)=(-1);
  g1(6,9)=1;
  g1(7,5)=(-1);
  g1(7,10)=1;
  g1(8,6)=(-1);
  g1(8,11)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,2744);
end
end
end
end
