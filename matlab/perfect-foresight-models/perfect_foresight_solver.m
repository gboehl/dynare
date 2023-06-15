function perfect_foresight_solver()
% Computes deterministic simulations
%
% INPUTS
%   None
%
% OUTPUTS
%   none
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 1996-2023 Dynare Team
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

global M_ options_ oo_ ys0_ ex0_

check_input_arguments(options_, M_, oo_);

periods = options_.periods;

if options_.debug
    model_static = str2func([M_.fname,'.static']);
    for ii=1:size(oo_.exo_simul,1)
        [residual(:,ii)] = model_static(oo_.steady_state, oo_.exo_simul(ii,:),M_.params);
    end
    problematic_periods=find(any(isinf(residual)) | any(isnan(residual)))-M_.maximum_endo_lag;
    if ~isempty(problematic_periods)
        period_string=num2str(problematic_periods(1));
        for ii=2:length(problematic_periods)
            period_string=[period_string, ', ', num2str(problematic_periods(ii))];
        end
        fprintf('\n\nWARNING: Value for the exogenous variable(s) in period(s) %s inconsistent with the static model.\n',period_string);
        fprintf('WARNING: Check for division by 0.\n')
    end
end

% Various sanity checks

if options_.no_homotopy && (options_.simul.homotopy_initial_step_size ~= 1 ...
                            || options_.simul.homotopy_max_completion_share ~= 1 ...
                            || options_.simul.homotopy_linearization_fallback ...
                            || options_.simul.homotopy_marginal_linearization_fallback ~= 0)
    error('perfect_foresight_solver: the no_homotopy option is incompatible with homotopy_initial_step_size, homotopy_max_completion_share, homotopy_linearization_fallback and homotopy_marginal_linearization_fallback options')
end

if options_.simul.homotopy_initial_step_size > 1 || options_.simul.homotopy_initial_step_size < 0
    error('perfect_foresight_solver: The value given to homotopy_initial_step_size option must be in the [0,1] interval')
end

if options_.simul.homotopy_min_step_size > 1 || options_.simul.homotopy_min_step_size < 0
    error('perfect_foresight_solver: The value given to homotopy_min_step_size option must be in the [0,1] interval')
end

if options_.simul.homotopy_max_completion_share > 1 || options_.simul.homotopy_max_completion_share < 0
    error('perfect_foresight_solver: The value given to homotopy_max_completion_share option must be in the [0,1] interval')
end

if options_.simul.homotopy_marginal_linearization_fallback > 1 || options_.simul.homotopy_marginal_linearization_fallback < 0
    error('perfect_foresight_solver: The value given to homotopy_marginal_linearization_fallback option must be in the [0,1] interval')
end

if options_.simul.homotopy_initial_step_size < options_.simul.homotopy_min_step_size
    error('perfect_foresight_solver: The value given to homotopy_initial_step_size option must be greater or equal to that given to homotopy_min_step_size option')
end

if options_.simul.homotopy_linearization_fallback && options_.simul.homotopy_marginal_linearization_fallback > 0
    error('perfect_foresight_solver: Options homotopy_linearization_fallback and homotopy_marginal_linearization_fallback cannot be used together')
end

if options_.simul.homotopy_max_completion_share < 1 && ~options_.simul.homotopy_linearization_fallback && options_.simul.homotopy_marginal_linearization_fallback == 0
    error('perfect_foresight_solver: Option homotopy_max_completion_share has a value less than 1, so you must also specify either homotopy_linearization_fallback or homotopy_marginal_linearization_fallback')
end

initperiods = 1:M_.maximum_lag;
simperiods = M_.maximum_lag+(1:periods);
lastperiods = M_.maximum_lag+periods+(1:M_.maximum_lead);

% Create base scenario for homotopy, which corresponds to the initial steady
% state (i.e. a known solution to the perfect foresight problem, assuming that
% oo_.steady_state/ys0_ effectively contains a steady state)
if isempty(ys0_)
    endobase = repmat(oo_.steady_state, 1,M_.maximum_lag+periods+M_.maximum_lead);
    exobase = repmat(oo_.exo_steady_state',M_.maximum_lag+periods+M_.maximum_lead,1);
else
    endobase = repmat(ys0_, 1, M_.maximum_lag+periods+M_.maximum_lead);
    exobase = repmat(ex0_', M_.maximum_lag+periods+M_.maximum_lead, 1);
end

% Determine whether the terminal condition is not a steady state (typically
% because steady was not called after endval)
if ~options_.simul.endval_steady && ~isempty(ys0_)
    terminal_condition_is_a_steady_state = true;
    for j = lastperiods
        endval_resid = evaluate_static_model(oo_.endo_simul(:,j), oo_.exo_simul(j,:), M_.params, M_, options_);
        if norm(endval_resid, 'Inf') > options_.simul.steady_tolf
            terminal_condition_is_a_steady_state = false;
            break
        end
    end
end

% Copy the paths for the exogenous and endogenous variables, as given by perfect_foresight_setup
exoorig = oo_.exo_simul;
endoorig = oo_.endo_simul;

current_share = 0;  % Share of shock successfully completed so far
step = min(options_.simul.homotopy_initial_step_size, options_.simul.homotopy_max_completion_share);
success_counter = 0;
iteration = 0;

function local_success = create_scenario(share)
    % For a given share, updates the exogenous path and also the initial and
    % terminal conditions for the endogenous path (but do not modify the initial
    % guess for endogenous)

    % Compute convex combination for the path of exogenous
    oo_.exo_simul = exoorig*share + exobase*(1-share);

    % Compute convex combination for the initial condition
    % In most cases, the initial condition is a steady state and this does nothing
    % This is for cases when the initial condition is out of equilibrium
    oo_.endo_simul(:, initperiods) = share*endoorig(:, initperiods)+(1-share)*endobase(:, initperiods);

    % If there is a permanent shock, ensure that the rescaled terminal condition is
    % a steady state (if the user asked for this recomputation, or if the original
    % terminal condition is a steady state)
    local_success = true;
    if options_.simul.endval_steady || (~isempty(ys0_) && terminal_condition_is_a_steady_state)

        % Set “local” options for steady state computation (after saving the global values)
        saved_steady_solve_algo = options_.solve_algo;
        options_.solve_algo = options_.simul.steady_solve_algo;
        saved_steady_maxit = options_.steady.maxit;
        options_.steady.maxit = options_.simul.steady_maxit;
        saved_steady_tolf = options_.solve_tolf;
        options_.solve_tolf = options_.simul.steady_tolf;
        saved_steady_tolx = options_.solve_tolx;
        options_.solve_tolx = options_.simul.steady_tolx;
        saved_steady_markowitz = options_.markowitz;
        options_.markowitz = options_.simul.steady_markowitz;

        saved_ss = oo_.endo_simul(:, lastperiods);
        % Effectively compute the terminal steady state
        for j = lastperiods
            % First use the terminal steady of the previous homotopy iteration as guess value (or the contents of the endval block if this is the first iteration)
            [oo_.endo_simul(:, j), ~, info] = evaluate_steady_state(oo_.endo_simul(:, j), oo_.exo_simul(j, :), M_, options_, true);
            if info(1)
                % If this fails, then try again using the initial steady state as guess value
                if isempty(ys0_)
                    guess_value = oo_.steady_state;
                else
                    guess_value = ys0_;
                end
                [oo_.endo_simul(:, j), ~, info] = evaluate_steady_state(guess_value, oo_.exo_simul(j, :), M_, options_, true);
                if info(1)
                    % If this fails again, give up and restore last periods in oo_.endo_simul
                    oo_.endo_simul(:, lastperiods) = saved_ss;
                    local_success = false;
                    break;
                end
            end
        end

        % The following is needed for the STEADY_STATE() operator to work properly
        oo_.steady_state = oo_.endo_simul(:, end);

        options_.solve_algo = saved_steady_solve_algo;
        options_.steady.maxit = saved_steady_maxit;
        options_.solve_tolf = saved_steady_tolf;
        options_.solve_tolx = saved_steady_tolx;
        options_.markowitz = saved_steady_markowitz;
    else
        % The terminal condition is not a steady state, compute a convex combination
        oo_.endo_simul(:, lastperiods) = share*endoorig(:, lastperiods)+(1-share)*endobase(:, lastperiods);
    end
end

while step > options_.simul.homotopy_min_step_size

    iteration = iteration+1;

    saved_endo_simul = oo_.endo_simul;

    new_share = current_share + step; % Try this share, and see if it succeeds

    if new_share > options_.simul.homotopy_max_completion_share
        new_share = options_.simul.homotopy_max_completion_share; % Don't go beyond target point
        step = new_share - current_share;
    end

    steady_success = create_scenario(new_share);

    if steady_success
        % At the first iteration, use the initial guess given by
        % perfect_foresight_setup or the user (but only if new_share=1, otherwise it
        % does not make much sense). Afterwards, until a converging iteration has been obtained,
        % use the rescaled terminal condition (or, if there is no lead, the base
        % scenario / initial steady state).
        if current_share == 0
            if iteration == 1 && new_share == 1
                oo_.endo_simul(:, simperiods) = endoorig(:, simperiods);
            elseif M_.maximum_lead > 0
                oo_.endo_simul(:, simperiods) = repmat(oo_.endo_simul(:, lastperiods(1)), 1, options_.periods);
            else
                oo_.endo_simul(:, simperiods) = endobase(:, simperiods);
            end
        end

        % Solve for the paths of the endogenous variables.
        [oo_.endo_simul, success, maxerror, solver_iter, per_block_status] = perfect_foresight_solver_core(M_, options_, oo_);
    else
        success = false;
        maxerror = NaN;
        solver_iter = [];
        per_block_status = [];
    end

    if options_.no_homotopy || (iteration == 1 && success && new_share == 1)
        % Skip homotopy
        if success
            current_share = new_share;
        end
        did_homotopy = false;
        break
    end

    if iteration == 1
        % First iteration failed, so we enter homotopy
        did_homotopy = true;

        if ~options_.noprint
            fprintf('\nEntering the homotopy method iterations...\n')
            fprintf('\nIter. \t | Share \t | Status \t | Max. residual\n')
            fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        end

        % Disable warnings if homotopy
        warning_old_state = warning;
        warning off all
        % Do not print anything
        oldverbositylevel = options_.verbosity;
        options_.verbosity = 0;
    end

    if success
        % Successful step
        if ~options_.noprint
            fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_share, 'succeeded', maxerror)
        end
        current_share = new_share;
        if current_share >= options_.simul.homotopy_max_completion_share
            break
        end
        success_counter = success_counter + 1;
        if options_.simul.homotopy_step_size_increase_success_count > 0 ...
           && success_counter >= options_.simul.homotopy_step_size_increase_success_count
            success_counter = 0;
            step = step * 2;
        end
    else
        oo_.endo_simul = saved_endo_simul;
        success_counter = 0;
        step = step / 2;
        if ~options_.noprint
            if ~steady_success
                fprintf('%i \t | %1.5f \t | %s\n', iteration, new_share, 'failed (in endval steady)')
            elseif isreal(maxerror)
                fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_share, 'failed', maxerror)
            else
                fprintf('%i \t | %1.5f \t | %s \t | %s\n', iteration, new_share, 'failed', 'Complex')
            end
        end
    end
end

if did_homotopy && ~options_.noprint
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')
end

%If simulated paths are complex, take real part and recompute the residuals to check whether this is actually a solution
if ~isreal(oo_.endo_simul(:)) % cannot happen with bytecode or the perfect_foresight_problem DLL
    real_simul = real(oo_.endo_simul);
    real_maxerror = recompute_maxerror(real_simul, oo_, M_, options_);
    if real_maxerror <= options_.dynatol.f
        oo_.endo_simul = real_simul;
        maxerror = real_maxerror;
    else
        current_share = 0;
        disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
    end
end

% Put solver status information in oo_, and do linearization if needed and requested
if current_share == 1
    if ~options_.noprint
        fprintf('Perfect foresight solution found.\n\n')
    end
    oo_.deterministic_simulation.status = true;
elseif options_.simul.homotopy_linearization_fallback && current_share > 0
    oo_.endo_simul = endobase + (oo_.endo_simul - endobase)/current_share;
    oo_.exo_simul = exoorig;
    if options_.simul.endval_steady
        % The following is needed for the STEADY_STATE() operator to work properly,
        % and thus must come before computing the maximum error.
        % This is not a true steady state, but it is the closest we can get to
        oo_.steady_state = oo_.endo_simul(:, end);
    end
    maxerror = recompute_maxerror(oo_.endo_simul, oo_, M_, options_);

    if ~options_.noprint
        fprintf('Perfect foresight solution found for %.1f%% of the shock, then extrapolation was performed using linearization\n\n', current_share*100)
    end
    oo_.deterministic_simulation.status = true;
    oo_.deterministic_simulation.homotopy_linearization = true;
elseif options_.simul.homotopy_marginal_linearization_fallback > 0 && current_share > options_.simul.homotopy_marginal_linearization_fallback
    saved_endo_simul = oo_.endo_simul;
    new_share = current_share - options_.simul.homotopy_marginal_linearization_fallback;
    new_success = create_scenario(new_share);
    if new_success
        [oo_.endo_simul, new_success] = perfect_foresight_solver_core(M_, options_, oo_);
    end
    if new_success
        oo_.endo_simul = saved_endo_simul + (saved_endo_simul - oo_.endo_simul)*(1-current_share)/options_.simul.homotopy_marginal_linearization_fallback;
        oo_.exo_simul = exoorig;
        if options_.simul.endval_steady
            % The following is needed for the STEADY_STATE() operator to work properly,
            % and thus must come before computing the maximum error.
            % This is not a true steady state, but it is the closest we can get to
            oo_.steady_state = oo_.endo_simul(:, end);
        end
        maxerror = recompute_maxerror(oo_.endo_simul, oo_, M_, options_);

        if ~options_.noprint
            fprintf('Perfect foresight solution found for %.1f%% of the shock, then extrapolation was performed using marginal linearization\n\n', current_share*100)
        end
        oo_.deterministic_simulation.homotopy_marginal_linearization = true;
    else
        fprintf('perfect_foresight_solver: marginal linearization failed, unable to find solution for %.1f%% of the shock. Try to modify the value of homotopy_marginal_linearization_fallback option\n\n', new_share*100)
    end
    oo_.deterministic_simulation.status = new_success;
else
    fprintf('Failed to solve perfect foresight model\n\n')
    oo_.deterministic_simulation.status = false;
end

oo_.deterministic_simulation.error = maxerror;
oo_.deterministic_simulation.homotopy_completion_share = current_share;
oo_.deterministic_simulation.homotopy_iterations = iteration;
if ~isempty(solver_iter)
    oo_.deterministic_simulation.iterations = solver_iter;
end
if ~isempty(per_block_status)
    oo_.deterministic_simulation.block = per_block_status;
end

% Must come after marginal linearization
if did_homotopy
    options_.verbosity = oldverbositylevel;
    warning(warning_old_state);
end

dyn2vec(M_, oo_, options_);

if isfield(oo_, 'initval_series') && ~isempty(oo_.initval_series)
    initial_period = oo_.initval_series.dates(1)+(M_.orig_maximum_lag-1);
elseif ~isdates(options_.initial_period) && isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    initial_period = options_.initial_period;
end

ts = dseries([transpose(oo_.endo_simul(1:M_.orig_endo_nbr,:)), oo_.exo_simul], initial_period, [M_.endo_names(1:M_.orig_endo_nbr); M_.exo_names]);

if isfield(oo_, 'initval_series') && ~isempty(oo_.initval_series)
    names = ts.name;
    ts = merge(oo_.initval_series{names{:}}, ts);
end

assignin('base', 'Simulated_time_series', ts);

if success
    oo_.gui.ran_perfect_foresight = true;
end

end

function maxerror = recompute_maxerror(endo_simul, oo_, M_, options_)
    % Computes ∞-norm of residuals for a given path of endogenous,
    % given the exogenous path and parameters in oo_ and M_
    if options_.bytecode
        residuals = bytecode('dynamic', 'evaluate', endo_simul, oo_.exo_simul, M_.params, oo_.steady_state, 1);
    else
        ny = size(endo_simul, 1);
        periods = size(endo_simul, 2) - M_.maximum_lag - M_.maximum_lead;
        if M_.maximum_lag > 0
            y0 = endo_simul(:, M_.maximum_lag);
        else
            y0 = NaN(ny, 1);
        end
        if M_.maximum_lead > 0
            yT = endo_simul(:, M_.maximum_lag+periods+1);
        else
            yT = NaN(ny, 1);
        end
        yy = endo_simul(:,M_.maximum_lag+(1:periods));
        residuals = perfect_foresight_problem(yy(:), y0, yT, oo_.exo_simul, M_.params, oo_.steady_state, periods, M_, options_);
    end
    maxerror = norm(vec(residuals), 'Inf');
end
