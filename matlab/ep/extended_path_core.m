function [y, info_convergence, endogenousvariablespaths] = extended_path_core(periods,endo_nbr,exo_nbr,positive_var_indx, ...
                                                  exo_simul,init,initial_conditions,...
                                                  steady_state, ...
                                                  debug,order,M_,pfm,algo,solve_algo,stack_solve_algo,...
                                                  olmmcp,options_,oo_,initialguess)

% Copyright © 2016-2023 Dynare Team
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

ep = options_.ep;

if init% Compute first order solution (Perturbation)...
    endo_simul = simult_(M_,options_,initial_conditions,oo_.dr,exo_simul(2:end,:),1);
else
    if nargin==19 && ~isempty(initialguess)
        % Note that the first column of initialguess should be equal to initial_conditions.
        endo_simul = initialguess;
    else
        endo_simul = [initial_conditions repmat(steady_state,1,periods+1)];
    end
end

oo_.endo_simul = endo_simul;

if debug
    save ep_test_1.mat endo_simul exo_simul
end

if options_.bytecode && order > 0
    error('Option order > 0 of extended_path command is not compatible with bytecode option.')
end
if options_.block && order > 0
    error('Option order > 0 of extended_path command is not compatible with block option.')
end

if order == 0
    options_.periods = periods;
    options_.block = pfm.block;
    oo_.endo_simul = endo_simul;
    oo_.exo_simul = exo_simul;
    oo_.steady_state = steady_state;
    options_.lmmcp = olmmcp;
    options_.solve_algo = solve_algo;
    options_.stack_solve_algo = stack_solve_algo;
    [endogenousvariablespaths, info_convergence] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
else
    switch(algo)
        case 0
            [flag, endogenousvariablespaths] = ...
            solve_stochastic_perfect_foresight_model(endo_simul, exo_simul, pfm, ep.stochastic.quadrature.nodes, ep.stochastic.order);
        case 1
            [flag, endogenousvariablespaths] = ...
            solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, options_, pfm, ep.stochastic.order);
    end
    info_convergence = ~flag;
end

if ~info_convergence && ~options_.no_homotopy
    [info_convergence, endogenousvariablespaths] = extended_path_homotopy(endo_simul, exo_simul, M_, options_, oo_, pfm, ep, order, algo, 2, debug);
end

if info_convergence
    y = endogenousvariablespaths(:,2);
else
    y = NaN(size(endo_nbr,1));
end
