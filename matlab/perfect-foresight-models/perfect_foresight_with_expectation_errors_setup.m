function perfect_foresight_with_expectation_errors_setup

% Copyright (C) 2021 Dynare Team
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

global M_ oo_ options_ ys0_

if ~isempty(ys0_)
    error('perfect_foresight_with_expectation_errors_setup: cannot be used in conjunction with endval')
end

if ~isempty(M_.endo_histval)
    error('perfect_foresight_with_expectation_errors_setup: cannot be used in conjunction with histval')
end
if ~isempty(oo_.initval_series)
    error('perfect_foresight_with_expectation_errors_setup: cannot be used in conjunction with histval_file/initval_file')
end

periods = options_.periods;

%% Read CSV file
if ~exist(options_.datafile, 'file')
    error(['perfect_foresight_with_expectation_errors_setup: cannot find ' options_.datafile ])
end
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
oo_.pfwee.terminal_info = NaN(M_.exo_nbr, periods); % 2nd dimension is informational time
oo_.pfwee.shocks_info = NaN(M_.exo_nbr, periods, periods); % 2nd dimension is real time, 3rd dimension is informational time
for i = 1:size(raw_csv, 2)
    exo_id = strmatch(exo_header_names{i}, M_.exo_names, 'exact');
    period_id = raw_csv(1, i);
    oo_.pfwee.terminal_info(exo_id, period_id) = raw_csv(2, i);
    % Ignore irrelevant periods when copying shocks information
    oo_.pfwee.shocks_info(exo_id, period_id:end, period_id) = raw_csv(2+period_id:end, i);
end

% Build initial paths for endos and exos (only initial conditions are set, the rest is NaN)
oo_.endo_simul = [repmat(oo_.steady_state, 1, M_.maximum_lag) NaN(M_.endo_nbr, periods+M_.maximum_lead)];
oo_.exo_simul = [repmat(oo_.exo_steady_state',M_.maximum_lag,1); NaN(periods+M_.maximum_lead,M_.exo_nbr)];
