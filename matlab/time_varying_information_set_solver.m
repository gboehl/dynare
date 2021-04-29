function time_varying_information_set_solver(csv_file, terminal_steady_state_as_guess_value)

global M_ oo_ options_ ys0_

if ~ischar(csv_file)
    error('A CSV file must be given as first argument')
end
if ~islogical(terminal_steady_state_as_guess_value)
    error('A boolean must be given as second argument')
end

if ~isempty(ys0_)
    error('Cannot be used in conjunction with endval')
end

periods = options_.periods;

%% Read CSV file
%% We can’t use readcell (only in MATLAB ≥ R2019a), so instead rely on csvread and manual hacks
% Read numeric data, skipping first row and first column
raw_csv = csvread(csv_file, 1, 1);
if size(raw_csv, 1)-2 ~= periods
    error(['The number of rows in ' csv_file ' does not match the periods setting'])
end
% Read first line (exogenous variable names)
fid = fopen(csv_file);
csv_first_line = fgetl(fid);
fclose(fid);
exo_header_names = strsplit(csv_first_line, ',');
exo_header_names = exo_header_names(2:end); % Remove first column
if numel(exo_header_names) ~= size(raw_csv, 2)
    error(['First line malformed in ' csv_file])
end

%% Create and fill structures containing information sets
terminal_info = NaN(M_.exo_nbr, periods); % 2nd dimension is informational time
shocks_info = NaN(M_.exo_nbr, periods, periods); % 2nd dimension is real time, 3rd dimension is informational time
for i = 1:size(raw_csv, 2)
    exo_id = strmatch(exo_header_names{i}, M_.exo_names, 'exact');
    period_id = raw_csv(1, i);
    terminal_info(exo_id, period_id) = raw_csv(2, i);
    % Ignore irrelevant periods when copying shocks information
    shocks_info(exo_id, period_id:end, period_id) = raw_csv(2+period_id:end, i);
end

if ~isempty(M_.endo_histval)
    error('histval unsupported')
end
if ~isempty(oo_.initval_series)
    error('histval_file/initval_file unsupported')
end

% Build initial paths for endos and exos (only initial conditions are set, the rest is NaN)
endo_simul = [repmat(oo_.steady_state, 1, M_.maximum_lag) NaN(M_.endo_nbr, periods+M_.maximum_lead)];
exo_simul = [repmat(oo_.exo_steady_state',M_.maximum_lag,1); NaN(periods+M_.maximum_lead,M_.exo_nbr)];

% Save initial steady state, for restoring it at the end
initial_steady_state = oo_.steady_state;
initial_exo_steady_state = oo_.exo_steady_state;

% Start main loop around informational periods
info_period = 1;
while info_period <= periods
    % Compute terminal steady state as anticipated
    oo_.exo_steady_state = terminal_info(:, info_period);
    steady_state_prev = oo_.steady_state;
    [oo_.steady_state,~,info] = evaluate_steady_state(steady_state_prev, M_, options_, oo_, true);

    options_.periods = periods - info_period + 1;
    oo_.endo_simul = endo_simul(:, info_period:end); % Take initial conditions + guess values from previous simulation
    if terminal_steady_state_as_guess_value
        % Overwrite guess value with terminal steady state
        oo_.endo_simul(:, M_.maximum_lag+(1:periods-info_period+1)) = repmat(oo_.steady_state, 1, periods-info_period+1);
    elseif info_period == 1
        % Use initial steady state as guess value for first simulation if not using terminal steady state
        oo_.endo_simul(:, M_.maximum_lag+(1:periods)) = repmat(initial_steady_state, 1, periods);
    end
    oo_.endo_simul(:, end-M_.maximum_lead+1:end) = repmat(oo_.steady_state, 1, M_.maximum_lead);
    oo_.exo_simul = exo_simul(info_period:end, :);
    oo_.exo_simul(M_.maximum_lag+(1:periods-info_period+1), :) = shocks_info(:, info_period:end, info_period)';
    oo_.exo_simul(end-M_.maximum_lead+1:end, :) = repmat(oo_.exo_steady_state, M_.maximum_lead, 1);

    perfect_foresight_solver;

    endo_simul(:, info_period:end) = oo_.endo_simul;
    exo_simul(info_period:end, :) = oo_.exo_simul;

    % Increment info_period (as much as possible, if information set does not change for some time)
    increment = 1;
    while info_period+increment <= periods && ...
          all(terminal_info(:, info_period) == terminal_info(:, info_period+increment)) && ...
          all(all(shocks_info(:, info_period+increment:end, info_period) == shocks_info(:, info_period+increment:end, info_period+increment)))
        increment = increment + 1;
    end
    info_period = info_period + increment;
end

options_.periods = periods;
oo_.endo_simul = endo_simul;
oo_.exo_simul = exo_simul;

oo_.steady_state = initial_steady_state;
oo_.exo_steady_state = initial_exo_steady_state;

end
