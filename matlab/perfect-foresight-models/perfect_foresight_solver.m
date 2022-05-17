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

% Copyright (C) 1996-2020 Dynare Team
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

global M_ options_ oo_

check_input_arguments(options_, M_, oo_);

if isempty(options_.scalv) || options_.scalv == 0
    options_.scalv = oo_.steady_state;
end

periods = options_.periods;

options_.scalv= 1;

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

initperiods = 1:M_.maximum_lag;
lastperiods = (M_.maximum_lag+periods+1):(M_.maximum_lag+periods+M_.maximum_lead);

oo_ = perfect_foresight_solver_core(M_,options_,oo_);

% If simulation failed try homotopy.
if ~oo_.deterministic_simulation.status && ~options_.no_homotopy

    if ~options_.noprint
        fprintf('\nSimulation of the perfect foresight model failed!')
        fprintf('Switching to a homotopy method...\n')
    end

    if ~M_.maximum_lag && M_.maximum_lead>0
        disp('Homotopy not implemented for purely forward models!')
        disp('Failed to solve the model!')
        disp('Return with empty oo_.endo_simul.')
        oo_.endo_simul = [];
        return
    end

    if ~M_.maximum_lead && M_.maximum_lag>0
        disp('Homotopy not implemented for purely backward models!')
        disp('Failed to solve the model!')
        disp('Return with empty oo_.endo_simul.')
        oo_.endo_simul = [];
        return
    end

    if ~M_.maximum_lead && ~M_.maximum_lag
        disp('Homotopy not implemented for purely static models!')
        disp('Failed to solve the model!')
        disp('Return with empty oo_.endo_simul.')
        oo_.endo_simul = [];
        return
    end

    % Disable warnings if homotopy
    warning_old_state = warning;
    warning off all
    % Do not print anything
    oldverbositylevel = options_.verbosity;
    options_.verbosity = 0;

    % Set initial paths for the endogenous and exogenous variables.
    endoinit = repmat(oo_.steady_state, 1,M_.maximum_lag+periods+M_.maximum_lead);
    exoinit = repmat(oo_.exo_steady_state',M_.maximum_lag+periods+M_.maximum_lead,1);

    % Copy the current paths for the exogenous and endogenous variables.
    exosim = oo_.exo_simul;
    endosim = oo_.endo_simul;

    current_weight = 0;    % Current weight of target point in convex combination.
    step = .5;             % Set default step size.
    success_counter = 0;
    iteration = 0;

    if ~options_.noprint
        fprintf('Iter. \t | Lambda \t | status \t | Max. residual\n')
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    end
    while (step > options_.dynatol.x)

        if ~isequal(step,1)
            options_.verbosity = 0;
        end

        iteration = iteration+1;
        new_weight = current_weight + step; % Try this weight, and see if it succeeds

        if new_weight >= 1
            new_weight = 1; % Don't go beyond target point
            step = new_weight - current_weight;
        end

        % Compute convex combination for exo path and initial/terminal endo conditions
        % But take care of not overwriting the computed part of oo_.endo_simul
        oo_.exo_simul = exosim*new_weight + exoinit*(1-new_weight);
        oo_.endo_simul(:,[initperiods, lastperiods]) = new_weight*endosim(:,[initperiods, lastperiods])+(1-new_weight)*endoinit(:,[initperiods, lastperiods]);

        % Detect Nans or complex numbers in the solution.
        path_with_nans = any(any(isnan(oo_.endo_simul)));
        path_with_cplx = any(any(~isreal(oo_.endo_simul)));

        if isequal(iteration, 1)
            % First iteration, same initial guess as in the first call to perfect_foresight_solver_core routine.
            oo_.endo_simul(:,M_.maximum_lag+1:end-M_.maximum_lead) = endoinit(:,1:periods);
        elseif path_with_nans || path_with_cplx
            % If solver failed with NaNs or complex number, use previous solution as an initial guess.
            oo_.endo_simul(:,M_.maximum_lag+1:end-M_.maximum_lead) = saved_endo_simul(:,1+M_.maximum_lag:end-M_.maximum_lead);
        end

        % Make a copy of the paths.
        saved_endo_simul = oo_.endo_simul;

        % Solve for the paths of the endogenous variables.
        [oo_,me] = perfect_foresight_solver_core(M_,options_,oo_);

        if oo_.deterministic_simulation.status == 1
            current_weight = new_weight;
            if current_weight >= 1
                if ~options_.noprint
                    fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'succeeded', me)
                end
                break
            end
            success_counter = success_counter + 1;
            if success_counter >= 3
                success_counter = 0;
                step = step * 2;
            end
            if ~options_.noprint
                fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'succeeded', me)
            end
        else
            % If solver failed, then go back.
            oo_.endo_simul = saved_endo_simul;
            success_counter = 0;
            step = step / 2;
            if ~options_.noprint
                if isreal(me)
                    fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'failed', me)
                else
                    fprintf('%i \t | %1.5f \t | %s \t | %s\n', iteration, new_weight, 'failed', 'Complex')
                end
            end
        end
    end
    if ~options_.noprint
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')
    end
    options_.verbosity = oldverbositylevel;
    warning(warning_old_state);
end

%If simulated paths are complex, take real part and recompute the residuals to check whether this is actually a solution
if ~isreal(oo_.endo_simul(:)) % cannot happen with bytecode or the perfect_foresight_problem DLL
    ny = size(oo_.endo_simul, 1);
    if M_.maximum_lag > 0
        y0 = real(oo_.endo_simul(:, M_.maximum_lag));
    else
        y0 = NaN(ny, 1);
    end
    if M_.maximum_lead > 0
        yT = real(oo_.endo_simul(:, M_.maximum_lag+periods+1));
    else
        yT = NaN(ny, 1);
    end
    if M_.maximum_lag~=0 && M_.maximum_lead~=0
        yy = real(oo_.endo_simul(:,M_.maximum_lag+(1:periods)));
        residuals = perfect_foresight_problem(yy(:), y0, yT, oo_.exo_simul, M_.params, oo_.steady_state, periods, M_, options_);
    else
        %The perfect_foresight_problem MEX only works on models with lags and leads
        i_cols = find(M_.lead_lag_incidence');
        residuals=NaN(ny,periods);
        yy=real(oo_.endo_simul);
        for it = (M_.maximum_lag+1):(M_.maximum_lag+periods)
            residuals(:,it) = feval([M_.fname '.dynamic'],yy(i_cols), oo_.exo_simul, M_.params, oo_.steady_state, it);
            i_cols = i_cols + ny;
        end
        residuals=residuals(:);
    end

    if max(abs(residuals))< options_.dynatol.f
        oo_.deterministic_simulation.status = 1;
        oo_.endo_simul=real(oo_.endo_simul);
    else
        oo_.deterministic_simulation.status = 0;
        disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
    end
end

if oo_.deterministic_simulation.status == 1
    if ~options_.noprint
        fprintf('Perfect foresight solution found.\n\n')
    end
else
    fprintf('Failed to solve perfect foresight model\n\n')
end

dyn2vec(M_, oo_, options_);

if ~isdates(options_.initial_period) && isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    initial_period = options_.initial_period;
end

ts = dseries(transpose(oo_.endo_simul), initial_period, M_.endo_names);
assignin('base', 'Simulated_time_series', ts);
if oo_.deterministic_simulation.status
    oo_.gui.ran_perfect_foresight = true;
end
