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

% Copyright (C) 2003-2020 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


% dseries
if isfield(options, 'series')
    series = options.series;
    dseries_ispresent = true;
else
    dseries_ispresent = false;
end

% file
datafile = '';
if isfield(options, 'filename')
    warning([caller, '_FILE: option FILENAME is deprecated, please use', ...
                 ' option DATAFILE'])
    if dseries_ispresent
        error([caller, '_FILE: you can''t use option FILENAME and option SERIES', ...
               ' at the same time'])
    end
    if isfield(options, 'datafile')
        error([caller, '_FILE: you can''t use option DATAFILE and option FILENAME', ...
               ' at the same time'])
    end
    datafile = options.filename;
end

if isfield(options, 'datafile')
    if dseries_ispresent
        error([caller, '_FILE: you can''t use option DATAFILE and option SERIES', ...
               ' at the same time'])
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
            error([caller, '_FILE: Can''t find datafile: ' basename '.{m,mat,xls,xlsx}']);
        end
    end

    fullname = [basename extension];
    series = dseries(fullname);
end

% checking that all variable are present
error_flag = false;
for i = 1:M.orig_endo_nbr
    if ~series.exist(M.endo_names{i})
        disp(sprintf('%s_FILE: endogenous variable %s is missing', ...
                     caller, M.endo_names{i}))
        error_flag = true;
    end
end

for i = 1:M.exo_nbr
    if ~series.exist(M.exo_names{i})
        disp(sprintf('%s_FILE: exogenous variable %s is missing', ...
                     caller, M.exo_names{i}))
        error_flag = true;
    end
end

for i = 1:M.exo_det_nbr
    if ~series.exist(M.exo_det_names{i})
        disp(sprintf('%s_FILE: exo_det variable %s is missing', ...
                     caller, M.exo_det_names{i}))
        error_flag = true;
    end
end

if error_flag
    error([caller, '_FILE: some variables are missing'])
end

if exist(sprintf('+%s/dynamic_set_auxiliary_series', M.fname), 'file')
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

if ~isfield(options, 'first_obs') || isempty(options.first_obs)
    if isfield(options, 'first_simulation_period')
        options.first_obs = options.first_simulation_period ...
            - M.orig_maximum_lag;
    end
elseif isfield(options, 'first_simulation_period')
    nobs = options.first_simulation_period - opions_.first_obs;
    if M.orig_maximum_lag ~= nobs
        error(sprintf(['HISTVALF: first_obs = %d and', ...
                       ' first_simulation_period = %d', ...
                       ' don''t provide for the number of' ...
                       ' lags in the model.'], ...
                      options.first_obs, ...
                      options.first_simulation_period))
    end
end

if isfield(options, 'first_obs')
    i = options.first_obs;
    if i < 1
        error([caller, '_FILE: the first requested period is before available', ...
               ' data.'])
    elseif i > nobs0
        error([caller, '_FILE: the first requested period is after available', ...
               ' data.'])
    end
    first_obs = periods(i);
    if nobs > 0
        last_obs = first_obs + nobs - 1;
        last_obs_ispresent = true;
    end
    first_obs_ispresent = true;
elseif isfield(options, 'firstobs')
    first_obs = options.firstobs; 
    if nobs > 0
        last_obs = first_obs + nobs - 1;
        last_obs_ispresent = true;
    end
    first_obs_ispresent = true;
end

if last_obs_ispresent
    if isfield(options, 'last_obs')
        i = options.last_obs;
        if i < 1
            error([caller, '_FILE: the last requested period is before available', ...
                   ' data.'])
        elseif i > nobs0
            error([caller, '_FILE: the last requested period is after available', ...
                   ' data.'])
        end
        if last_obs ~= periods(i) 
            error([caller, '_FILE: FIST_OBS, LAST_OBS and NOBS contain', ...
                   ' inconsistent information. Use only two of these', ...
                   ' options.'])
        end    
    elseif isfield(options, 'lastobs')
        if last_obs ~= options.lastobs 
            error([caller, '_FILE: FIST_OBS, LAST_OBS and NOBS contain', ...
                   ' inconsistent information. Use only two of these', ...
                   ' options.'])
        end    
    end
elseif isfield(options, 'last_obs')
    i = options.last_obs;
    if i < 1
        error([caller, '_FILE: the last requested period is before available', ...
               ' data.'])
    elseif i > nobs0
        error([caller, '_FILE: the last requested period is after available', ...
               ' data.'])
    end
    last_obs = periods(i);
    if nobs > 0
        first_obs = last_obs - nobs + 1;
        first_obs_ispresent = true;
    end
    last_obs_ispresent = true;
elseif isfield(options, 'lastobs')
    last_obs = options.lastobs; 
    if nobs > 0
        first_obs = last_obs - nobs + 1;
        first_obs_ispresent = true;
    end
    last_obs_ispresent = true;
end

if ~first_obs_ispresent
    first_obs = periods(1);
end

if ~last_obs_ispresent
    if nobs > 0
        last_obs = first_obs + nobs - 1;
    else
        last_obs = periods(end);
    end
end

if first_obs < series.init
    error([caller, '_FILE: the first requested period is before available', ...
                    ' data.'])
elseif last_obs > series.last
    error([caller, '_FILE: the last requested period is after available', ...
                    ' data.'])
else
    series = series(first_obs:last_obs);
end

