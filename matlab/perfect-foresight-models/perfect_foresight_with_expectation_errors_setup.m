function oo_=perfect_foresight_with_expectation_errors_setup(M_, options_, oo_)
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

if ~isempty(M_.endo_histval)
    error('perfect_foresight_with_expectation_errors_setup: cannot be used in conjunction with histval')
end
if ~isempty(oo_.initval_series)
    error('perfect_foresight_with_expectation_errors_setup: cannot be used in conjunction with histval_file/initval_file')
end

periods = options_.periods;

%% Initialize informational structures
oo_.pfwee.terminal_info = NaN(M_.exo_nbr, periods); % 2nd dimension is informational time
oo_.pfwee.shocks_info = NaN(M_.exo_nbr, periods, periods); % 2nd dimension is real time, 3rd dimension is informational time

if exist(options_.datafile, 'file')
    if ~isempty(M_.det_shocks) || ~isempty(M_.learnt_shocks) || ~isempty(oo_.initial_steady_state) || ~isempty(M_.learnt_endval)
        warning('perfect_foresight_with_expectation_errors_setup: since you passed the datafile option, the contents of shocks and endval blocks will be ignored')
    end
    %% Read CSV file
    %% We can’t use readcell (only in MATLAB ≥ R2019a), so instead rely on csvread and manual hacks
    % Read numeric data, skipping first row and first column
    raw_csv = csvread(options_.datafile, 1, 1);
    if size(raw_csv, 1)-2 ~= periods
        error(['perfect_foresight_with_expectation_errors_setup: the number of rows in ' options_.datafile ' does not match the periods setting'])
    end
    % Read first line (exogenous variable names)
    fid = fopen(options_.datafile);
    csv_first_line = fgetl(fid);
    fclose(fid);
    exo_header_names = strsplit(csv_first_line, ',');
    exo_header_names = exo_header_names(2:end); % Remove first column
    if numel(exo_header_names) ~= size(raw_csv, 2)
        error(['perfect_foresight_with_expectation_errors_setup: first line malformed in ' options_.datafile])
    end

    %% Create and fill structures containing information sets
    for i = 1:size(raw_csv, 2)
        exo_id = strmatch(exo_header_names{i}, M_.exo_names, 'exact');
        period_id = raw_csv(1, i);
        % Ignore irrelevant periods when copying shocks information
        oo_.pfwee.shocks_info(exo_id, period_id:end, period_id) = raw_csv(1+period_id:end-1, i);
        oo_.pfwee.terminal_info(exo_id, period_id) = raw_csv(end, i);
    end
else
    %% No datafile option given, use the contents of shocks and endval blocks
    if isempty(M_.learnt_shocks) && isempty(M_.learnt_endval)
        warning('perfect_foresight_with_expectation_errors_setup: there is no shocks(learnt_in=...) or endval(learnt_in=...) block, and you did not pass the datafile option, so there is no point in using this command')
    end

    %% Initialize information set at period 1 using “bare” shocks and endval blocks (or initval if there is no endval)
    oo_.pfwee.terminal_info(:, 1) = oo_.exo_steady_state;
    oo_.pfwee.shocks_info(:, :, 1) = repmat(oo_.exo_steady_state, 1, periods);
    for i = 1:length(M_.det_shocks)
        prds = M_.det_shocks(i).periods;
        exo_id = M_.det_shocks(i).exo_id;
        v = M_.det_shocks(i).value;
        if ~M_.det_shocks(i).exo_det
            switch M_.det_shocks(i).type
                case 'level'
                    oo_.pfwee.shocks_info(exo_id, prds, 1) = v;
                case 'multiply_steady_state'
                    oo_.pfwee.shocks_info(exo_id, prds, 1) = oo_.exo_steady_state(exo_id) * v;
                case 'multiply_initial_steady_state'
                    if isempty(oo_.initial_exo_steady_state)
                        error('Option relative_to_initval of mshocks block cannot be used without an endval block')
                    end
                    oo_.pfwee.shocks_info(exo_id, prds, 1) = oo_.initial_exo_steady_state(exo_id) * v;
            end
        end
    end

    %% Construct information sets for subsequent informational periods
    for p = 2:periods
        oo_.pfwee.terminal_info(:, p) = oo_.pfwee.terminal_info(:, p-1);
        oo_.pfwee.shocks_info(:, :, p) = oo_.pfwee.shocks_info(:, :, p-1);
        if ~isempty(M_.learnt_endval)
            idx = find([M_.learnt_endval.learnt_in] == p);
            for i = 1:length(idx)
                j = idx(i);
                exo_id = M_.learnt_endval(j).exo_id;
                switch M_.learnt_endval(j).type
                    case 'level'
                        oo_.pfwee.terminal_info(exo_id, p) = M_.learnt_endval(j).value;
                    case 'add'
                        oo_.pfwee.terminal_info(exo_id, p) = oo_.pfwee.terminal_info(exo_id, p-1) + M_.learnt_endval(j).value;
                    case 'multiply'
                        oo_.pfwee.terminal_info(exo_id, p) = oo_.pfwee.terminal_info(exo_id, p-1) * M_.learnt_endval(j).value;
                    otherwise
                        error('Unknown type in M_.learnt_endval')
                end
                oo_.pfwee.shocks_info(exo_id, p:end, p) = oo_.pfwee.terminal_info(exo_id, p);
            end
        end
        if ~isempty(M_.learnt_shocks)
            idx = find([M_.learnt_shocks.learnt_in] == p);
            for i = 1:length(idx)
                j = idx(i);
                exo_id = M_.learnt_shocks(j).exo_id;
                prds = M_.learnt_shocks(j).periods;
                switch M_.learnt_shocks(j).type
                    case 'level'
                        oo_.pfwee.shocks_info(exo_id, prds, p) = M_.learnt_shocks(j).value;
                    case 'add'
                        oo_.pfwee.shocks_info(exo_id, prds, p) = oo_.pfwee.shocks_info(exo_id, prds, p-1) + M_.learnt_shocks(j).value;
                    case 'multiply'
                        oo_.pfwee.shocks_info(exo_id, prds, p) = oo_.pfwee.shocks_info(exo_id, prds, p-1) .* M_.learnt_shocks(j).value;
                    case 'multiply_steady_state'
                        oo_.pfwee.shocks_info(exo_id, prds, p) = oo_.pfwee.terminal_info(exo_id, p) * M_.learnt_shocks(j).value;
                    otherwise
                        error('Unknown type in M_.learnt_shocks')
                end
            end
        end
    end
end

% Build initial paths for endos and exos (only initial conditions are set, the rest is NaN)
if isempty(oo_.initial_steady_state)
    oo_.endo_simul = repmat(oo_.steady_state, 1, M_.maximum_lag+periods+M_.maximum_lead);
else
    oo_.endo_simul = [repmat(oo_.initial_steady_state, 1, M_.maximum_lag) repmat(oo_.steady_state, 1, periods+M_.maximum_lead)];
end
if isempty(oo_.initial_exo_steady_state)
    oo_.exo_simul = repmat(oo_.exo_steady_state', M_.maximum_lag+periods+M_.maximum_lead, 1);
else
    oo_.exo_simul = [repmat(oo_.initial_exo_steady_state', M_.maximum_lag, 1); repmat(oo_.exo_steady_state', periods+M_.maximum_lead, 1)];
end
