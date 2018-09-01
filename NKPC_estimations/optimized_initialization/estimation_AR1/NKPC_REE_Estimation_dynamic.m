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
y_forecast__ = params(11)+params(11)*params(13)+params(13)^2*y(1);
pinf_forecast__ = params(12)+params(12)*params(14)+params(14)^2*y(2);
lhs =y(6);
rhs =y_forecast__-1/params(5)*(y(8)-pinf_forecast__)+y(9);
residual(1)= lhs-rhs;
lhs =y(7);
rhs =pinf_forecast__*beta__+y(6)*params(4)+y(10);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =params(10)*y(3)+(1-params(10))*(y(7)*params(6)+y(6)*params(7))+x(it_, 3);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =params(8)*y(4)+x(it_, 1);
residual(4)= lhs-rhs;
lhs =y(10);
rhs =params(9)*y(5)+x(it_, 2);
residual(5)= lhs-rhs;
lhs =y(11);
rhs =y(6)+params(1);
residual(6)= lhs-rhs;
lhs =y(12);
rhs =y(7)+params(2);
residual(7)= lhs-rhs;
lhs =y(13);
rhs =y(8)+params(3);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 16);

  %
  % Jacobian matrix
  %

T21 = params(14)^2;
  g1(1,1)=(-(params(13)^2));
  g1(1,6)=1;
  g1(1,2)=1/params(5)*(-T21);
  g1(1,8)=1/params(5);
  g1(1,9)=(-1);
  g1(2,6)=(-params(4));
  g1(2,2)=(-(T21*beta__));
  g1(2,7)=1;
  g1(2,10)=(-1);
  g1(3,6)=(-((1-params(10))*params(7)));
  g1(3,7)=(-((1-params(10))*params(6)));
  g1(3,3)=(-params(10));
  g1(3,8)=1;
  g1(3,16)=(-1);
  g1(4,4)=(-params(8));
  g1(4,9)=1;
  g1(4,14)=(-1);
  g1(5,5)=(-params(9));
  g1(5,10)=1;
  g1(5,15)=(-1);
  g1(6,6)=(-1);
  g1(6,11)=1;
  g1(7,7)=(-1);
  g1(7,12)=1;
  g1(8,8)=(-1);
  g1(8,13)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,256);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,4096);
end
end
end
end
