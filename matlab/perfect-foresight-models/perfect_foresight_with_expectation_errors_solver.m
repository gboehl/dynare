function oo_=perfect_foresight_with_expectation_errors_solver(M_, options_, oo_)
% INPUTS
%   M_                  [structure] describing the model
%   options_            [structure] describing the options
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   oo_                 [structure] storing the results

% Copyright © 2021-2023 Dynare Team
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

% Same for periods (it will be modified before calling perfect_foresight_solver if constants_simulation_length option is false)
periods = options_.periods;

% Retrieve initial paths built by pfwee_setup
% (the versions in oo_ will be truncated before calling perfect_foresight_solver)
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;

% Enforce the endval steady option in pf_solver
if M_.maximum_lead > 0
    options_.simul.endval_steady = true;
end

% Start main loop around informational periods
info_period = 1;
increment = 0;
while info_period <= periods
    if ~options_.noprint
        fprintf('perfect_foresight_with_expectations_errors_solver: computing solution for information available at period %d\n', info_period)
    end

    % Compute terminal steady state as anticipated
    oo_.exo_steady_state = oo_.pfwee.terminal_info(:, info_period);
    % oo_.steady_state will be updated by pf_solver (since endval_steady=true)

    if options_.pfwee.constant_simulation_length && increment > 0
        % Use previous terminal steady state as guess value for: simulation periods that don’t yet have an initial guess (i.e. are NaNs at this point); and also for the terminal steady state
        endo_simul = [ endo_simul repmat(oo_.steady_state, 1, increment)];
        exo_simul = [ exo_simul; NaN(increment, M_.exo_nbr)];
    end

    oo_.endo_simul = endo_simul(:, info_period:end); % Take initial conditions + guess values from previous simulation
    if options_.pfwee.constant_simulation_length
        sim_length = periods;
    else
        sim_length = periods - info_period + 1;
    end

    oo_.exo_simul = exo_simul(info_period:end, :);
    oo_.exo_simul(M_.maximum_lag+(1:periods-info_period+1), :) = oo_.pfwee.shocks_info(:, info_period:end, info_period)';
    oo_.exo_simul(M_.maximum_lag+periods-info_period+2:end, :) = repmat(oo_.exo_steady_state', sim_length+M_.maximum_lead-(periods-info_period+1), 1);

    options_.periods = sim_length;

    if info_period > 1 && homotopy_completion_share < 1 && options_.simul.homotopy_marginal_linearization_fallback > 0
        marginal_linearization_previous_raw_sims.sim1.endo_simul = oo_.deterministic_simulation.sim1.endo_simul(:, info_period:end);
        marginal_linearization_previous_raw_sims.sim1.exo_simul = oo_.deterministic_simulation.sim1.exo_simul(info_period:end, :);
        marginal_linearization_previous_raw_sims.sim1.homotopy_completion_share = oo_.deterministic_simulation.sim1.homotopy_completion_share;
        marginal_linearization_previous_raw_sims.sim2.endo_simul = oo_.deterministic_simulation.sim2.endo_simul(:, info_period:end);
        marginal_linearization_previous_raw_sims.sim2.exo_simul = oo_.deterministic_simulation.sim2.exo_simul(info_period:end, :);
        marginal_linearization_previous_raw_sims.sim2.homotopy_completion_share = oo_.deterministic_simulation.sim2.homotopy_completion_share;
    else
        marginal_linearization_previous_raw_sims = [];
    end

    oo_= perfect_foresight_solver(M_, options_, oo_, true, marginal_linearization_previous_raw_sims);

    if ~oo_.deterministic_simulation.status
        error('perfect_foresight_with_expectation_errors_solver: failed to compute solution for information available at period %d\n', info_period)
    end

    if info_period == 1
        homotopy_completion_share = oo_.deterministic_simulation.homotopy_completion_share;
        options_.simul.homotopy_max_completion_share = homotopy_completion_share;
    elseif oo_.deterministic_simulation.homotopy_completion_share ~= homotopy_completion_share
        error('perfect_foresight_solver_with_expectation_errors: could not find a solution for information available at period %d with the same homotopy completion share as period 1\n', info_period)
    end

    endo_simul(:, info_period:end) = oo_.endo_simul;
    exo_simul(info_period:end, :) = oo_.exo_simul;

    % Increment info_period (as much as possible, if information set does not change for some time)
    increment = 1;
    while info_period+increment <= periods && ...
          all(oo_.pfwee.terminal_info(:, info_period) == oo_.pfwee.terminal_info(:, info_period+increment)) && ...
          all(all(oo_.pfwee.shocks_info(:, info_period+increment:end, info_period) == oo_.pfwee.shocks_info(:, info_period+increment:end, info_period+increment)))
        increment = increment + 1;
    end
    info_period = info_period + increment;
end

% Set final paths
oo_.endo_simul = endo_simul;
oo_.exo_simul = exo_simul;