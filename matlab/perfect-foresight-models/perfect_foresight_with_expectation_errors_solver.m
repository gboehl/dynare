function perfect_foresight_with_expectation_errors_solver

% Copyright © 2021-2022 Dynare Team
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

global M_ oo_ options_

% Save initial steady state, for restoring it at the end
initial_steady_state = oo_.steady_state;
initial_exo_steady_state = oo_.exo_steady_state;

% Same for periods (it will be modified before calling perfect_foresight_solver if constants_simulation_length option is false)
periods = options_.periods;

% Retrieve initial paths built by pfwee_setup
% (the versions in oo_ will be truncated before calling perfect_foresight_solver)
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;

% Start main loop around informational periods
info_period = 1;
increment = 0;
while info_period <= periods
    % Compute terminal steady state as anticipated
    oo_.exo_steady_state = oo_.pfwee.terminal_info(:, info_period);
    steady_state_prev = oo_.steady_state;
    [oo_.steady_state,~,info] = evaluate_steady_state(steady_state_prev, M_, options_, oo_, true);

    if options_.pfwee.constant_simulation_length && increment > 0
        endo_simul = [ endo_simul NaN(M_.endo_nbr, increment)];
        exo_simul = [ exo_simul; NaN(increment, M_.exo_nbr)];
    end

    oo_.endo_simul = endo_simul(:, info_period:end); % Take initial conditions + guess values from previous simulation
    if options_.pfwee.constant_simulation_length
        sim_length = periods;
    else
        sim_length = periods - info_period + 1;
    end
    if options_.pfwee.terminal_steady_state_as_guess_value
        % Overwrite guess value with terminal steady state
        oo_.endo_simul(:, M_.maximum_lag+(1:sim_length)) = repmat(oo_.steady_state, 1, sim_length);
    elseif info_period == 1
        % Use initial steady state as guess value for first simulation if not using terminal steady state
        oo_.endo_simul(:, M_.maximum_lag+(1:periods)) = repmat(initial_steady_state, 1, periods);
    elseif options_.pfwee.constant_simulation_length && increment > M_.maximum_lead
        % Use initial steady state as guess value for simulation periods that don’t yet have an initial guess (i.e. are NaNs)
        oo_.endo_simul(:, M_.maximum_lag+periods-(0:increment-M_.maximum_lead-1)) = repmat(initial_steady_state, 1, increment-M_.maximum_lead);
    end
    oo_.endo_simul(:, end-M_.maximum_lead+1:end) = repmat(oo_.steady_state, 1, M_.maximum_lead);
    oo_.exo_simul = exo_simul(info_period:end, :);
    oo_.exo_simul(M_.maximum_lag+(1:periods-info_period+1), :) = oo_.pfwee.shocks_info(:, info_period:end, info_period)';
    oo_.exo_simul(M_.maximum_lag+periods-info_period+2:end) = repmat(oo_.exo_steady_state, sim_length+M_.maximum_lead-(periods-info_period+1), 1);

    options_.periods = sim_length;

    perfect_foresight_solver;

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

% Restore some values
oo_.steady_state = initial_steady_state;
oo_.exo_steady_state = initial_exo_steady_state;
options_.periods = periods;
