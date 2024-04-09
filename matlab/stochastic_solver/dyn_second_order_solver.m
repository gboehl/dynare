function dr = dyn_second_order_solver(g1, g2, dr, M_, threads_BC)

%@info:
%! @deftypefn {Function File} {@var{dr} =} dyn_second_order_solver (@var{g1}, @var{g2}, @var{dr}, @var{M_}, @var{threads_BC})
%! @anchor{dyn_second_order_solver}
%! @sp 1
%! Computes the second order reduced form of the DSGE model, for details please refer to
%! * Juillard and Kamenik (2004): Solving Stochastic Dynamic Equilibrium Models: A k-Order Perturbation Approach
%! * Kamenik (2005) - Solving SDGE Models: A New Algorithm for the Sylvester Equation
%! Note that this function makes use of the fact that Dynare internally transforms the model
%! so that there is only one lead and one lag on endogenous variables and, in the case of a stochastic model,
%! no leads/lags on exogenous variables. See the manual for more details.
%  Auxiliary variables
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item g1
%! Sparse matrix containing the Jacobian of the dynamic model
%! @item g2
%! Sparse matrix containing the Hessian of the dynamic model
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item M_
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item threads_BC
%! Integer controlling number of threads in sparse_hessian_times_B_kronecker_C
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @end table
%! @end deftypefn
%@eod:

% Copyright © 2001-2024 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

dr.ghxx = [];
dr.ghuu = [];
dr.ghxu = [];
dr.ghs2 = [];

% Indices in the DR-reordered variables
klag = M_.nstatic+(1:M_.nspred);
[~, kcurr] = find(M_.lead_lag_incidence(2, dr.order_var));
klead = M_.nstatic+M_.npred+(1:M_.nsfwrd);

% Indices in the g1 columns that map to DR-order
ilag = dr.order_var(klag);
icurr = M_.endo_nbr + dr.order_var(kcurr);
ilead = 2*M_.endo_nbr + dr.order_var(klead);

kk1 = [ilag; icurr; ilead;
       3*M_.endo_nbr+(1:M_.exo_nbr+M_.exo_det_nbr)'];
nk = size(g1, 2);
kk2 = reshape(1:nk^2,nk,nk);

%% ghxx
A = zeros(M_.endo_nbr,M_.endo_nbr);
A(:,kcurr) = g1(:,icurr);
A(:,klag) = A(:,klag) + g1(:,ilead)*dr.ghx(klead,:);
B = zeros(M_.endo_nbr,M_.endo_nbr);
B(:,M_.nstatic+M_.npred+1:end) = g1(:, ilead);
C = dr.ghx(klag,:);
zx = [eye(M_.nspred);
      dr.ghx(kcurr,:);
      dr.ghx(klead,:)*dr.ghx(klag,:);
      zeros(M_.exo_nbr,M_.nspred);
      zeros(M_.exo_det_nbr,M_.nspred)];
zu = [zeros(M_.nspred,M_.exo_nbr);
      dr.ghu(kcurr,:);
      dr.ghx(klead,:)*dr.ghu(klag,:);
      eye(M_.exo_nbr);
      zeros(M_.exo_det_nbr,M_.exo_nbr)];
rhs = sparse_hessian_times_B_kronecker_C(g2(:,kk2(kk1,kk1)), zx, threads_BC); % g2: reordering to DR order
rhs = -rhs;
dr.ghxx = gensylv(2,A,B,C,rhs);


%% ghxu
%rhs
rhs = sparse_hessian_times_B_kronecker_C(g2(:,kk2(kk1,kk1)), zx, zu, threads_BC); % g2: reordering to DR order
abcOut = A_times_B_kronecker_C(dr.ghxx, dr.ghx(klag,:), dr.ghu(klag,:));
rhs = -rhs-B*abcOut;
%lhs
dr.ghxu = A\rhs;

%% ghuu
%rhs
rhs = sparse_hessian_times_B_kronecker_C(g2(:,kk2(kk1,kk1)), zu, threads_BC); % g2: reordering to DR order
B1 = A_times_B_kronecker_C(B*dr.ghxx,dr.ghu(klag,:));
rhs = -rhs-B1;
%lhs
dr.ghuu = A\rhs;

%% ghs2
% derivatives of F with respect to forward variables
O1 = zeros(M_.endo_nbr,M_.nstatic);
O2 = zeros(M_.endo_nbr,M_.nfwrd);
LHS = zeros(M_.endo_nbr,M_.endo_nbr);
LHS(:,kcurr) = g1(:,icurr);
RHS = zeros(M_.endo_nbr,M_.exo_nbr^2);
E = eye(M_.endo_nbr);
B1 = sparse_hessian_times_B_kronecker_C(g2(:,kk2(ilead,ilead)), dr.ghu(klead,:), threads_BC); % g2: focus only on forward variables and reorder to DR order
RHS = RHS + g1(:,ilead)*dr.ghuu(klead,:)+B1;
% LHS
LHS = LHS + g1(:,ilead)*(E(klead,:)+[O1(klead,:) dr.ghx(klead,:) O2(klead,:)]);
RHS = RHS*M_.Sigma_e(:);
dr.fuu = RHS;
RHS = -RHS;
dr.ghs2 = LHS\RHS;
