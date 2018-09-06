function [residual, g1, g2, g3] = NKPC_REE_Simulation_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(5, 1);
beta__ = 0.99;
lhs =y(4);
rhs =y(9)-1/params(5)*(y(6)-y(10))+y(7);
residual(1)= lhs-rhs;
lhs =y(5);
rhs =y(10)*beta__+y(4)*params(4)+y(8);
residual(2)= lhs-rhs;
lhs =y(6);
rhs =params(10)*y(1)+(1-params(10))*(y(5)*params(6)+y(4)*params(7))+x(it_, 3);
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(8)*y(2)+x(it_, 1);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(9)*y(3)+x(it_, 2);
residual(5)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(5, 13);

  %
  % Jacobian matrix
  %

  g1(1,4)=1;
  g1(1,9)=(-1);
  g1(1,10)=(-(1/params(5)));
  g1(1,6)=1/params(5);
  g1(1,7)=(-1);
  g1(2,4)=(-params(4));
  g1(2,5)=1;
  g1(2,10)=(-beta__);
  g1(2,8)=(-1);
  g1(3,4)=(-((1-params(10))*params(7)));
  g1(3,5)=(-((1-params(10))*params(6)));
  g1(3,1)=(-params(10));
  g1(3,6)=1;
  g1(3,13)=(-1);
  g1(4,2)=(-params(8));
  g1(4,7)=1;
  g1(4,11)=(-1);
  g1(5,3)=(-params(9));
  g1(5,8)=1;
  g1(5,12)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],5,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],5,2197);
end
end
end
end
