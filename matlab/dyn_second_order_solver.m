function dr = dyn_second_order_solver(jacobia,hessian_mat,dr,M,threads_BC)

%@info:
%! @deftypefn {Function File} {@var{dr} =} dyn_second_order_solver (@var{jacobia},@var{hessian_mat},@var{dr},@var{M_},@var{threads_BC})
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
%! @item jacobia
%! Matrix containing the Jacobian of the model
%! @item hessian_mat
%! Matrix containing the second order derivatives of the model
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

% Copyright (C) 2001-2020 Dynare Team
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

k1 = nonzeros(M.lead_lag_incidence(:,dr.order_var)');
kk1 = [k1; length(k1)+(1:M.exo_nbr+M.exo_det_nbr)'];
nk = size(kk1,1);
kk2 = reshape(1:nk^2,nk,nk);
ic = [ M.nstatic+(1:M.nspred) ]';

klag  = M.lead_lag_incidence(1,dr.order_var); %columns are in DR order
kcurr = M.lead_lag_incidence(2,dr.order_var); %columns are in DR order
klead = M.lead_lag_incidence(3,dr.order_var); %columns are in DR order

%% ghxx
A = zeros(M.endo_nbr,M.endo_nbr);
A(:,kcurr~=0) = jacobia(:,nonzeros(kcurr));
A(:,ic) = A(:,ic) + jacobia(:,nonzeros(klead))*dr.ghx(klead~=0,:);
B = zeros(M.endo_nbr,M.endo_nbr);
B(:,M.nstatic+M.npred+1:end) = jacobia(:,nonzeros(klead));
C = dr.ghx(ic,:);
zx = [eye(length(ic));
      dr.ghx(kcurr~=0,:);
      dr.ghx(klead~=0,:)*dr.ghx(ic,:);
      zeros(M.exo_nbr,length(ic));
      zeros(M.exo_det_nbr,length(ic))];
zu = [zeros(length(ic),M.exo_nbr);
      dr.ghu(kcurr~=0,:);
      dr.ghx(klead~=0,:)*dr.ghu(ic,:);
      eye(M.exo_nbr);
      zeros(M.exo_det_nbr,M.exo_nbr)];
rhs = sparse_hessian_times_B_kronecker_C(hessian_mat(:,kk2(kk1,kk1)),zx,threads_BC); %hessian_mat: reordering to DR order
rhs = -rhs;
dr.ghxx = gensylv(2,A,B,C,rhs);


%% ghxu
%rhs
rhs = sparse_hessian_times_B_kronecker_C(hessian_mat(:,kk2(kk1,kk1)),zx,zu,threads_BC); %hessian_mat: reordering to DR order
abcOut = A_times_B_kronecker_C(dr.ghxx, dr.ghx(ic,:), dr.ghu(ic,:));
rhs = -rhs-B*abcOut;
%lhs
dr.ghxu = A\rhs;

%% ghuu
%rhs
rhs = sparse_hessian_times_B_kronecker_C(hessian_mat(:,kk2(kk1,kk1)),zu,threads_BC); %hessian_mat: reordering to DR order
B1 = A_times_B_kronecker_C(B*dr.ghxx,dr.ghu(ic,:));
rhs = -rhs-B1;
%lhs
dr.ghuu = A\rhs;

%% ghs2
% derivatives of F with respect to forward variables
O1 = zeros(M.endo_nbr,M.nstatic);
O2 = zeros(M.endo_nbr,M.nfwrd);
LHS = zeros(M.endo_nbr,M.endo_nbr);
LHS(:,kcurr~=0) = jacobia(:,nonzeros(kcurr));
RHS = zeros(M.endo_nbr,M.exo_nbr^2);
E = eye(M.endo_nbr);
B1 = sparse_hessian_times_B_kronecker_C(hessian_mat(:,kk2(nonzeros(klead),nonzeros(klead))), dr.ghu(klead~=0,:),threads_BC); %hessian_mat:focus only on forward variables and reorder to DR order
RHS = RHS + jacobia(:,nonzeros(klead))*dr.ghuu(klead~=0,:)+B1;
% LHS
LHS = LHS + jacobia(:,nonzeros(klead))*(E(klead~=0,:)+[O1(klead~=0,:) dr.ghx(klead~=0,:) O2(klead~=0,:)]);
RHS = RHS*M.Sigma_e(:);
dr.fuu = RHS;
RHS = -RHS;
dr.ghs2 = LHS\RHS;
