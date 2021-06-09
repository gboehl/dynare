function series = histvalf_initvalf(caller, M, options) 
% function initvalf(M)
%
% handles options for histvalf_initvalf() and initvalf()
%
% INPUTS
%    caller:           string, name of calling function
%    M:                model structure
%    options:          options specific to initivalf
%
% OUTPUTS
%    series:           dseries containing selected data from a file or a dseries
%

% Copyright (C) 2003-2021 Dynare Team
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


% dseries
if isfield(options, 'series')
    series = evalin('base', options.series);
    dseries_ispresent = true;
else
    dseries_ispresent = false;
end

% file
datafile = '';
if isfield(options, 'filename')
    warning('%s_FILE: option FILENAME is deprecated, please use option DATAFILE', caller)
    if dseries_ispresent
        error('%s_FILE: you can''t use option FILENAME and option SERIES at the same time', caller)
    end
    if isfield(options, 'datafile')
        error('%s_FILE: you can''t use option DATAFILE and option FILENAME at the same time', caller)
    end
    datafile = options.filename;
end

if isfield(options, 'datafile')
    if dseries_ispresent
        error('%s_FILE: you can''t use option DATAFILE and option SERIES at the same time', caller)
    end
    datafile = options.datafile;
end

if datafile
    [directory,basename,extension] = fileparts(datafile);
    % Auto-detect extension if not provided
    if isempty(extension)
        if exist([basename '.m'],'file')
            extension = '.m';
        elseif exist([basename '.mat'],'file')
            extension = '.mat';
        elseif exist([basename '.xls'],'file')
            extension = '.xls';
        elseif exist([basename '.xlsx'],'file')
            extension = '.xlsx';
        else
            error('%s_FILE: Can''t find datafile: %s.{m,mat,xls,xlsx}', caller, basename);
        end
    end
    fullname = [basename extension];
    series = dseries(fullname);
end

% checking that all variable are present
error_flag = false;
for i = 1:M.orig_endo_nbr
    if ~series.exist(M.endo_names{i})
        dprintf('%s_FILE: endogenous variable %s is missing', caller, M.endo_names{i})
        error_flag = true;
    end
end

for i = 1:M.exo_nbr
    if ~series.exist(M.exo_names{i})
        dprintf('%s_FILE: exogenous variable %s is missing', caller, M.exo_names{i})
        error_flag = true;
    end
end

for i = 1:M.exo_det_nbr
    if ~series.exist(M.exo_det_names{i})
        dprintf('%s_FILE: exo_det variable %s is missing', caller, M.exo_det_names{i})
        error_flag = true;
    end
end

if error_flag
    error('%s_FILE: some variables are missing', caller)
end

if exist(sprintf('+%s/dynamic_set_auxiliary_series.m', M.fname), 'file')
    series = feval(sprintf('%s.dynamic_set_auxiliary_series', M.fname), series, M.params);
end

% selecting observations
if isfield(options, 'nobs')
    nobs = options.nobs;
else
    nobs = 0;
end

periods = series.dates;
nobs0 = series.nobs;

first_obs_ispresent = false;
last_obs_ispresent = false;

first_obs = periods(1);
if isfield(options, 'first_obs') && ~isempty(options.first_obs)
    if options.first_obs < 1
        error('%s_FILE: first_obs must be a positive number', caller)
    elseif options.first_obs > nobs0
        error('%s_FILE: first_obs = %d is larger than the number of observations in the data file (%d)', ...
                      caller, options.first_obs, nobs0)
    elseif isfield(options, 'first_simulation_period')
        if  options.first_obs == options.first_simulation_period - M.orig_maximum_lag
            first_obs = periods(options.first_obs);
        else
            error('%s_FILE: first_obs = %d and first_simulation_period = %d have values inconsistent with a maximum lag of %d periods', ...
                          caller, options.first_obs, options.first_simulation_period, M.orig_maximum_lag)
        end
    elseif isfield(options, 'firstsimulationperiod')
        if  periods(options.first_obs) == options.firstsimulationperiod - M.orig_maximum_lag
            first_obs = periods(options.first_obs);
        else
            error('%s_FILE: first_obs = %d and first_simulation_period = %s have values inconsistent with a maximum lag of %d periods', ...
                          caller, options.first_obs, options.firstsimulationperiod, M.orig_maximum_lag)
        end
    else
        first_obs = periods(options.first_obs);
    end
    first_obs_ispresent = true;
end

if isfield(options, 'firstobs') && ~isempty(options.firstobs)
    if isfield(options, 'first_simulation_period')
        if  options.firstobs == periods(options.first_simulation_period) - M.orig_maximum_lag
            first_obs = options.firstobs;
        else
            error('%s_FILE: first_obs = %s and first_simulation_period = %d have values inconsistent with a maximum lag of %d periods', ...
                          caller, options.firstobs, options.first_simulation_period, M.orig_maximum_lag)
        end
    elseif isfield(options, 'firstsimulationperiod')
        if  options.firstobs == options.firstsimulationperiod - M.orig_maximum_lag
            first_obs = options.firstobs;
        else
            error('%s_FILE: firstobs = %s and first_simulation_period = %s have values inconsistent with a maximum lag of %d periods', ...
                          caller, options.firstobs, options.firstsimulationperiod, M.orig_maximum_lag)
        end
    else
        first_obs = options.firstobs;
    end
    first_obs_ispresent = true;
end

if ~first_obs_ispresent
    if isfield(options, 'first_simulation_period')
        if options.first_simulation_period < M.orig_maximum_lag
            error('%s_FILE: first_simulation_period = %d must be larger than the maximum lag (%d)', ...
                          caller, options.first_simulation_period, M.orig_maximum_lag)
        elseif options.first_simulation_period > nobs0
            error('%s_FILE: first_simulations_period = %d is larger than the number of observations in the data file (%d)', ...
                          caller, options.first_obs, nobs0)
        else
            first_obs = periods(options.first_simulation_period) - M.orig_maximum_lag;
        end
        first_obs_ispresent = true;
    elseif isfield(options, 'firstsimulationperiod')
        first_obs = options.firstsimulationperiod - M.orig_maximum_lag;
        first_obs_ispresent = true;
    end
end

if isfield(options, 'last_obs')
    if options.last_obs > nobs0
        error('%s_FILE: last_obs = %d is larger than the number of observations in the dataset (%d)', ...
                      caller, options.last_obs, nobs0)
    elseif first_obs_ispresent
        if nobs > 0 && (periods(options.last_obs) ~= first_obs + nobs - 1)
            error('%s_FILE: FIST_OBS, LAST_OBS and NOBS contain inconsistent information. Use only two of these options.', caller)
        else
            last_obs = periods(options.last_obs);
        end
    else
        last_obs = periods(options.last_obs);
        if nobs > 0
            first_obs = last_obs - nobs + 1;
        else
            first_obs = periods(1);
        end
    end
elseif isfield(options, 'lastobs')
    if options.lastobs > series.last
        error('%s_FILE: last_obs = %s is larger than the number of observations in the dataset (%s)', ...
                      caller, options.lastobs, series.last)
    elseif first_obs_ispresent
        if nobs > 0 && (options.lastobs ~= first_obs + nobs - 1)
            error('%s_FILE: FIST_OBS, LAST_OBS and NOBS contain inconsistent information. Use only two of these options.', caller)
        else
            last_obs = options.lastobs;
        end
    else
        last_obs = options.last_obs;
        if nobs > 0
            first_obs = last_obs - nobs + 1;
        else
            first_obs = periods(1);
        end
    end
elseif nobs > 0
    last_obs = first_obs + nobs - 1;
else
    last_obs = series.last;
end

series = series(first_obs:last_obs);
