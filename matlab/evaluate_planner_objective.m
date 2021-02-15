function planner_objective_value = evaluate_planner_objective(M_,options_,oo_)

% OUTPUT
% Returns a vector containing first order or second-order approximations of 
% - the unconditional expectation of the planner's objective function
% - the conditional expectation of the planner's objective function starting from the non-stochastic steady state and allowing for future shocks
% depending on the value of options_.order.

% ALGORITHM
% Welfare verifies
% W(y_{t-1}, u_t, sigma) = U(h(y_{t-1}, u_t, sigma)) + beta E_t W(g(y_{t-1}, u_t, sigma), u_t, sigma)
% where
% - W is the welfare function
% - U is the utility function
% - y_{t-1} is the vector of state variables
% - u_t is the vector of exogenous shocks scaled with sigma i.e. u_t = sigma e_t where e_t is the vector of exogenous shocks
% - sigma is the perturbation parameter
% - h is the policy function, providing controls x_t in function of information at time t i.e. (y_{t-1}, u_t, sigma)
% - g is the transition function, providing next-period state variables in function of information at time t i.e. (y_{t-1}, u_t, sigma)
% - beta is the planner's discount factor
% - E_t is the expectation operator given information at time t i.e. (y_{t-1}, u_t, sigma)

% The unconditional expectation of the planner's objective function verifies
% E(W) = E(U)/(1-beta)
% The conditional expectation of the planner's objective function given (y_{t-1}, u_t, sigma) coincides with the welfare function delineated above.

% A first-order approximation of the utility function around the non-stochastic steady state (y_{t-1}, u_t, sigma) = (y, 0, 0) is
% U(h(y_{t-1}, u_t, sigma)) = Ubar + U_x ( h_y yhat_{t-1} + h_u u_t )
% Taking the unconditional expectation yields E(U) = Ubar and E(W) = Ubar/(1-beta)
% As for conditional welfare, a first-order approximation leads to
% W = Wbar + W_y yhat_{t-1} + W_u u_t
% The approximated conditional expectation of the planner's objective function taking at the non-stochastic steady-state and allowing for future shocks thus verifies
% W (y, 0, 1) = Wbar

% Similarly, taking the unconditional expectation of a second-order approximation of utility around the non-stochastic steady state yields a second-order approximation of unconditional welfare
% E(W) = (1 - beta)^{-1} ( Ubar + U_x h_y E(yhat) + 0.5 ( (U_xx h_y^2 + U_x h_yy) E(yhat^2) + (U_xx h_u^2 + U_x h_uu) E(u^2) + U_x h_ss )
% where E(yhat), E(yhat^2) and E(u^2) can be derived from oo_.mean and oo_.var

% As for conditional welfare, the second-order approximation of welfare around the non-stochastic steady state leads to
% W(y_{t-1}, u_t, sigma) = Wbar + W_y yhat_{t-1} + W_u u_t + W_yu yhat_{t-1} ⊗ u_t + 0.5 ( W_yy yhat_{t-1}^2 + W_uu u_t^2 + W_ss )
% The derivatives of W taken at the non-stochastic steady state can be computed as in Kamenik and Juillard (2004) "Solving Stochastic Dynamic Equilibrium Models: A k-Order Perturbation Approach".
% The approximated conditional expectation of the planner's objective function starting from the non-stochastic steady-state and allowing for future shocks thus verifies
% W(y,0,1) = Wbar + 0.5*Wss

% INPUTS
%   M_:        (structure) model description
%   options_:  (structure) options
%   oo_:       (structure) output results
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2007-2021 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

dr = oo_.dr;

exo_nbr = M_.exo_nbr;
nstatic = M_.nstatic;
nspred = M_.nspred;
beta = get_optimal_policy_discount_factor(M_.params, M_.param_names);

ys = oo_.dr.ys;
planner_objective_value = zeros(2,1);
if options_.order == 1
    [U] = feval([M_.fname '.objective.static'],ys,zeros(1,exo_nbr), M_.params);
    planner_objective_value(1) = U/(1-beta);
    planner_objective_value(2) = U/(1-beta);  
elseif options_.order == 2
    [U,Uy,Uyy] = feval([M_.fname '.objective.static'],ys,zeros(1,exo_nbr), M_.params);

    Gy = dr.ghx(nstatic+(1:nspred),:);
    Gu = dr.ghu(nstatic+(1:nspred),:);
    Gyy = dr.ghxx(nstatic+(1:nspred),:);
    Gyu = dr.ghxu(nstatic+(1:nspred),:);
    Guu = dr.ghuu(nstatic+(1:nspred),:);
    Gss = dr.ghs2(nstatic+(1:nspred),:);

    gy(dr.order_var,:) = dr.ghx;
    gu(dr.order_var,:) = dr.ghu;
    gyy(dr.order_var,:) = dr.ghxx;
    gyu(dr.order_var,:) = dr.ghxu;
    guu(dr.order_var,:) = dr.ghuu;
    gss(dr.order_var,:) = dr.ghs2;

    Uyy = full(Uyy);

    Uyygygy = A_times_B_kronecker_C(Uyy,gy,gy);
    Uyygugu = A_times_B_kronecker_C(Uyy,gu,gu);

    %% Unconditional welfare

    old_noprint = options_.noprint;

    if ~old_noprint
        options_.noprint = 1;
    end
    var_list = M_.endo_names(dr.order_var(nstatic+(1:nspred)));
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list); %get decision rules and moments
    if ~old_noprint
        options_.noprint = 0;
    end

    oo_.mean(isnan(oo_.mean)) = options_.huge_number;
    oo_.var(isnan(oo_.var)) = options_.huge_number;

    Ey = oo_.mean;
    Eyhat = Ey - ys(dr.order_var(nstatic+(1:nspred)));

    var_corr = Eyhat*Eyhat';
    Eyhatyhat = oo_.var(:) + var_corr(:);
    Euu = M_.Sigma_e(:);

    EU = U + Uy*gy*Eyhat + 0.5*((Uyygygy + Uy*gyy)*Eyhatyhat + (Uyygugu + Uy*guu)*Euu + Uy*gss);
    EW = EU/(1-beta);

    %% Conditional welfare starting from the non-stochastic steady-state

    Wbar = U/(1-beta);
    Wy = Uy*gy/(eye(nspred)-beta*Gy);

    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;
    end
    %solve Lyapunuv equation Wyy=gy'*Uyy*gy+Uy*gyy+beta*Wy*Gyy+beta*Gy'Wyy*Gy
    Wyy = reshape(lyapunov_symm(sqrt(beta)*Gy',reshape(Uyygygy + Uy*gyy + beta*Wy*Gyy,nspred,nspred),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold, 3, options_.debug),1,nspred*nspred);
    Wyygugu = A_times_B_kronecker_C(Wyy,Gu,Gu);
    Wuu = Uyygugu + Uy*guu + beta*(Wyygugu + Wy*Guu);
    Wss = (Uy*gss + beta*(Wy*Gss + Wuu*M_.Sigma_e(:)))/(1-beta);
    W = Wbar + 0.5*Wss;

    planner_objective_value(1) = EW;
    planner_objective_value(2) = W;
else
    %Order k code will go here!
    fprintf('\nevaluate_planner_objective: order>2 not yet supported\n')
    planner_objective_value(1) = NaN;
    planner_objective_value(2) = NaN;
    return
end
if ~options_.noprint
    if options_.ramsey_policy
        fprintf('\nApproximated value of unconditional welfare:  %10.8f\n', planner_objective_value(1))
        fprintf('\nApproximated value of conditional welfare:  %10.8f\n', planner_objective_value(2))
    elseif options_.discretionary_policy
        fprintf('\nApproximated value of unconditional welfare with discretionary policy: %10.8f\n\n', EW)
    end
end
