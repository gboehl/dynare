function planner_objective_value = evaluate_planner_objective(M_,options_,oo_)
% function planner_objective_value = evaluate_planner_objective(M_,options_,oo_)
% INPUTS
%   M_:        (structure) model description
%   options_:  (structure) options
%   oo_:       (structure) output results
% OUTPUT
%  planner_objective_value (double)
%
%Returns a vector containing first order or second-order approximations of
% - the unconditional expectation of the planner's objective function
% - the conditional expectation of the planner's objective function starting from the non-stochastic steady state and allowing for future shocks
% depending on the value of options_.order.
%
% SPECIAL REQUIREMENTS
%   none

% ALGORITHM
% Welfare satifies
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

% The unconditional expectation of the planner's objective function satisfies
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

% In the discretionary case, the model is assumed to be linear and the utility is assumed to be linear-quadratic. This changes 2 aspects of the results delinated above:
% 1) the second-order derivatives of the policy and transition functions h and g are zero.
% 2) the unconditional expectation of states coincides with its steady-state, which entails E(yhat) = 0
% Therefore, the unconditional welfare can now be approximated as
% E(W) = (1 - beta)^{-1} ( Ubar + 0.5 ( U_xx h_y^2 E(yhat^2) + U_xx h_u^2 E(u^2) )
% As for the conditional welfare, the second-order formula above is still valid, but the derivatives of W no longer contain any second-order derivatives of the policy and transition functions h and g.

% In the deterministic case, resorting to approximations for welfare is no longer required as it is possible to simulate the model given initial conditions for pre-determined variables and terminal conditions for forward-looking variables, whether these initial and terminal conditions are explicitly or implicitly specified. Assuming that the number of simulated periods is high enough for the new steady-state to be reached, the new unconditional welfare is thus the last period's welfare. As for the conditional welfare, it can be derived using backward recursions on the equation W = U + beta*W(+1) starting from the final unconditional steady-state welfare.

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

dr = oo_.dr;
if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

exo_nbr = M_.exo_nbr;
nstatic = M_.nstatic;
nspred = M_.nspred;
beta = get_optimal_policy_discount_factor(M_.params, M_.param_names);

planner_objective_value = zeros(2,1);
if options_.ramsey_policy
    if oo_.gui.ran_perfect_foresight
        T = size(oo_.endo_simul,2);
        [U_term] = feval([M_.fname '.objective.static'],oo_.endo_simul(:,T-M_.maximum_lead),oo_.exo_simul(T-M_.maximum_lead,:), M_.params);
        EW = U_term/(1-beta);
        W = EW;
        for t=T-M_.maximum_lead:-1:1+M_.maximum_lag
            [U] = feval([M_.fname '.objective.static'],oo_.endo_simul(:,t),oo_.exo_simul(t,:), M_.params);
            W = U + beta*W;
        end
        planner_objective_value(1) = EW;
        planner_objective_value(2) = W;
    else
        ys = oo_.dr.ys;
        if options_.order == 1 || M_.hessian_eq_zero
            [U] = feval([M_.fname '.objective.static'],ys,zeros(1,exo_nbr), M_.params);
            planner_objective_value(1) = U/(1-beta);
            planner_objective_value(2) = U/(1-beta);
        elseif options_.order == 2 && ~M_.hessian_eq_zero
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
            if options_.pruning
                fprintf('evaluate_planner_objective: pruning option is not supported and will be ignored\n')
            end
            oo_=disp_th_moments(dr,var_list,M_,options_,oo_);
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
            fprintf('\nevaluate_planner_objective: order>2 unconditional welfare calculation not yet supported\n')
            planner_objective_value(1) = k_order_welfare(dr, M_, options_);
            planner_objective_value(2) = NaN;
            return
        end
    end
elseif options_.discretionary_policy
    ys = oo_.dr.ys;
    
    [U,Uy,Uyy] = feval([M_.fname '.objective.static'],ys,zeros(1,exo_nbr), M_.params);
    
    Gy = dr.ghx(nstatic+(1:nspred),:);
    Gu = dr.ghu(nstatic+(1:nspred),:);
    gy(dr.order_var,:) = dr.ghx;
    gu(dr.order_var,:) = dr.ghu;
    
    Uyy = full(Uyy);
    
    Uyygygy = A_times_B_kronecker_C(Uyy,gy,gy);
    Uyygugu = A_times_B_kronecker_C(Uyy,gu,gu);
    
    %% Unconditional welfare
    
    old_noprint = options_.noprint;
    
    if ~old_noprint
        options_.noprint = 1;
    end
    var_list = M_.endo_names(dr.order_var(nstatic+(1:nspred)));
    oo_=disp_th_moments(dr,var_list,M_,options_,oo_);
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
    
    EU = U + Uy*gy*Eyhat + 0.5*(Uyygygy*Eyhatyhat + Uyygugu*Euu);
    EW = EU/(1-beta);
    
    
    %% Conditional welfare starting from the non-stochastic steady-state
    
    Wbar = U/(1-beta);
    Wy = Uy*gy/(eye(nspred)-beta*Gy);
    
    %solve Lyapunuv equation Wyy=gy'*Uyy*gy+beta*Gy'Wyy*Gy
    Wyy = reshape(lyapunov_symm(sqrt(beta)*Gy',reshape(Uyygygy,nspred,nspred),options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold, 3, options_.debug),1,nspred*nspred);
    Wyygugu = A_times_B_kronecker_C(Wyy,Gu,Gu);
    Wuu = Uyygugu + beta*Wyygugu;
    Wss = beta*Wuu*M_.Sigma_e(:)/(1-beta);
    W = Wbar + 0.5*Wss;
    
    planner_objective_value(1) = EW;
    planner_objective_value(2) = W;
end
if ~options_.noprint
    if options_.ramsey_policy
        if oo_.gui.ran_perfect_foresight
            fprintf('\nSimulated value of unconditional welfare:  %10.8f\n', planner_objective_value(1))
            fprintf('\nSimulated value of conditional welfare:  %10.8f\n', planner_objective_value(2))
        else
            fprintf('\nApproximated value of unconditional welfare:  %10.8f\n', planner_objective_value(1))
            fprintf('\nApproximated value of conditional welfare:  %10.8f\n', planner_objective_value(2))
        end
    elseif options_.discretionary_policy
        fprintf('\nApproximated value of unconditional welfare with discretionary policy:  %10.8f\n', planner_objective_value(1))
        fprintf('\nApproximated value of conditional welfare with discretionary policy:  %10.8f\n', planner_objective_value(2))
    end
end
