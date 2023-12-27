function [dataset_, dataset_info, newdatainterface] = makedataset(options_, initialconditions, gsa_flag)
%[dataset_, dataset_info, newdatainterface] = makedataset(options_, initialconditions, gsa_flag)
% Initialize a dataset as a dseries object.
% INPUTS
% ======
%
%     options_              [struct]    Structure of options built by Dynare's preprocessor.
%     initialconditions     [double]    number of lags for VAR and DSGE_VAR
%     gsa_flag              [integer]   1: GSA, 0: other
%
% OUTPUTS
% =======
%
%     dataset_       [dseries]  The dataset.
%     dataset_info   [struct]   Various informations about the dataset (descriptive statistics and missing observations).
%
% EXAMPLE
% =======
%
%     [dataset_, dataset_info] = makedataset(options_) ;
%
%
% See also dynare_estimation_init

% Copyright © 2014-2023 Dynare Team
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


if nargin<3
    gsa_flag = 0;
end

if nargin<2 || isempty(initialconditions)
    % If a the sample is to be used for the estimation of a VAR or DSGE-VAR model
    % the second argument must be a strictly positive integer (the number of lags).
    initialconditions = 0;
end

if isempty(options_.datafile) && isempty(options_.dataset.file) && isempty(options_.dataset.series)
    if gsa_flag
        dataset_ = dseries();
        dataset_info = struct('missing', struct('state', 0, 'aindex', [], 'vindex', [], 'number_of_observations', NaN, 'no_more_missing_observations', NaN), ...
                             'descriptive', struct('mean', [], 'covariance', [], 'correlation', [], 'autocovariance', []));
        newdatainterface=0;
        return
    else
        error('makedataset: datafile option is missing!')
    end
end

if isempty(options_.datafile) && ~isempty(options_.dataset.file)
    datafile = options_.dataset.file;
    newdatainterface = 1;
elseif isempty(options_.datafile) && ~isempty(options_.dataset.series)
    try
        dseriesobjectforuserdataset = evalin('base', options_.dataset.series);
    catch
        error('makedataset: %s is unknown!', options_.dataset.series)
    end
    if ~isdseries(dseriesobjectforuserdataset)
        error('makedataset: %s has to be a dseries object!', options_.dataset.series)
    end
    datafile = [];
    newdatainterface = 1;
elseif ~isempty(options_.datafile) && isempty(options_.dataset.file)
    datafile = options_.datafile;
    newdatainterface = 0;
elseif ~isempty(options_.datafile) && ~isempty(options_.dataset.file)
    error('makedataset: You cannot simultaneously use the data command and the datafile option (in the estimation command)!')
else
    error('makedataset: You have to specify the datafile!')
end

% Check extension.
if ~isempty(datafile)
    allowed_extensions = {'m','mat','csv','xls','xlsx'};
    datafile_extension = get_file_extension(datafile);
    if isempty(datafile_extension)
        available_extensions = {}; j = 1;
        [datafilepath, datafilename] = fileparts(datafile);
        if isempty(datafilepath)
            datafilepath = '.';
        end
        dircontent = dir(datafilepath);
        for i=1:length(allowed_extensions)
            if ~isempty(strmatch([datafilename '.' allowed_extensions{i}],{dircontent.name},'exact'))
                available_extensions(j) = {allowed_extensions{i}};
                j = j+1;
            end
        end
        if isempty(available_extensions)
            error('makedataset: I can''t find a datafile (with allowed extension m, mat, csv, xls or xlsx)!')
        end
        if length(available_extensions)>1
            error(sprintf(['makedataset: You did not specify an extension for the datafile, but more than one candidate ' ...
                           'is available in the designated folder!\nPlease, add an extension to the datafile ' ...
                           '(m, mat, csv, xls or xlsx are permitted extensions).']));
        end
        datafile = [datafile '.' available_extensions{1}];
    end
end

% Load the data in a dseries object.
if ~isempty(datafile)
    if ~(newdatainterface==0 && ((length(datafile)>2 && strcmp(datafile(end-1:end),'.m')) || (length(datafile)>4 && strcmp(datafile(end-3:end),'.mat'))))
        dataset_ = dseries(datafile);
    else
        if length(datafile)>2 && strcmp(datafile(end-1:end),'.m')
            % Load an m file with the old interface.
            dataset_ = load_m_file_data_legacy(datafile, options_.varobs);
        elseif length(datafile)>4 && strcmp(datafile(end-3:end),'.mat')
            % Load a mat file with the old interface.
            dataset_ = load_mat_file_data_legacy(datafile, options_.varobs);
        end
    end
else
    dataset_ = dseriesobjectforuserdataset;
    clear('dseriesobjectforuserdataset');
end

if size(unique(dataset_.name),1)~=size(dataset_.name,1)
    error('makedataset: the data set must not contain two variables with the same name and must not contain empty/non-named columns.')
end

% Select a subset of the variables.
dataset_ = dataset_{options_.varobs{:}};

% Apply log function if needed.
if options_.loglinear && ~options_.logdata
    dataset_ = dataset_.log();
end

% Test if an initial period (different from its default value) is explicitely defined in the datafile.
if isequal(dataset_.init, dates(1,1))
    dataset_default_initial_period = 1;
else
    dataset_default_initial_period = 0;
end

%  Test if an initial period (different from its default value) is explicitely defined in the mod file with the set_time command.
if ~isdates(options_.initial_period) && isnan(options_.initial_period)
    set_time_default_initial_period = 1;
else
    set_time_default_initial_period = 0;
end

if ~set_time_default_initial_period && dataset_default_initial_period
    % Overwrite the initial period in dataset (it was set to default).
    % Note that the updates of freq and time members are auto-magically
    % done by dseries::subsasgn overloaded method.
    dataset_.init = options_.initial_period;
end

if set_time_default_initial_period && ~dataset_default_initial_period
    % Overwrite the global initial period defined by set_time (it was set to default).
    options_.initial_period = dataset_.init;
end

if ~set_time_default_initial_period && ~dataset_default_initial_period
    % Check if dataset.init and options_.initial_period are identical.
    if options_.initial_period<dataset_.init
        error('makedataset: The date as defined by the set_time command is not consistent with the initial period in the database!')
    end
end

% Set firstobs, lastobs and nobs
if newdatainterface
    if isempty(options_.dataset.firstobs)
        % first_obs option was not used in the data command.
        firstobs = dataset_.init;
    else
        firstobs = options_.dataset.firstobs;
    end
    if isnan(options_.dataset.nobs)
        % nobs option was not used in the data command.
        if isempty(options_.dataset.lastobs)
            % last_obs option was not used in the data command.
            nobs = dataset_.nobs;
            lastobs = dataset_.dates(end);
        else
            lastobs = options_.dataset.lastobs;
            nobs = lastobs-firstobs+1;
        end
    else
        nobs = options_.dataset.nobs;
        if isempty(options_.dataset.lastobs)
            % last_obs option was not used in the data command.
            lastobs = firstobs+(nobs-1);
        else
            % last_obs and nobs were used in the data command. Check that they are consistent (with firstobs).
            if ~isequal(lastobs,firstobs+(nobs-1))
                error('makedataset: Options last_obs (%s), first_obs (%s) and nobs (%s) are not consistent!',char(lastobs),char(firstobs),num2str(nobs));
            end
        end
    end
else
    if isnan(options_.first_obs)
        firstobs = dataset_.init;
    else
        firstobs = dataset_.dates(options_.first_obs);
    end
    if isnan(options_.nobs)
        lastobs = dataset_.dates(end);
        nobs = lastobs-firstobs+1;
    else
        nobs = options_.nobs;
        lastobs = firstobs+(nobs-1);
    end
end

% Add initial conditions if needed
FIRSTOBS = firstobs-initialconditions;

% Check that firstobs belongs to dataset_.dates
if firstobs<dataset_.init
    error('makedataset: first_obs (%s) cannot be less than the first date in the dataset (%s)!',char(firstobs),char(dataset_.init))
end

% Check that FIRSTOBS belongs to dataset_.dates
if initialconditions && FIRSTOBS<dataset_.init
    error('makedataset: first_obs (%s) - %i cannot be less than the first date in the dataset (%s)!\nReduce the number of lags in the VAR model or increase the value of first_obs\nto at least first_obs=%i.', char(firstobs), initialconditions, char(dataset_.init),initialconditions+1);
end

% Check that lastobs belongs to dataset_.dates...
if newdatainterface
    if lastobs>dataset_.dates(end)
        error('makedataset: last_obs (%s) cannot be greater than the last date in the dataset (%s)!',char(lastobs),char(dataset_.dates(end)))
    end
else
    % ...  or check that nobs is smaller than the number of observations in dataset_.
    if FIRSTOBS>dataset_.dates(1)
        if FIRSTOBS+nobs-1>dataset_.dates(end)
            error('makedataset: Given first_obs=%u and %u total observations in the dataset, the current nobs of %s must not be greater than %s!', options_.first_obs, dataset_.nobs, num2str(nobs), num2str(dataset_.nobs-find(dataset_.dates==FIRSTOBS)+1))
        end
    else
        if nobs>dataset_.nobs
            error('makedataset: nobs (%s) cannot be greater than the last date in the dataset (%s)!', num2str(nobs), num2str(dataset_.nobs))
        end
    end
end

% Select a subsample.
dataset_ = dataset_(FIRSTOBS:lastobs);

% Initialize dataset_info structure.
dataset_info = struct('missing', struct('state', NaN, 'aindex', [], 'vindex', [], 'number_of_observations', NaN, 'no_more_missing_observations', NaN), ...
                     'descriptive', struct('mean', [], 'covariance', [], 'correlation', [], 'autocovariance', []));

% Fill dataset_info.missing if some observations are missing
dataset_info.missing.state = isanynan(dataset_.data);
if dataset_info.missing.state
    [dataset_info.missing.aindex, dataset_info.missing.number_of_observations, dataset_info.missing.no_more_missing_observations, dataset_info.missing.vindex] = ...
        describe_missing_data(dataset_.data);
else
    dataset_info.missing.aindex = num2cell(transpose(repmat(1:dataset_.vobs,dataset_.nobs,1)),1);
    dataset_info.missing.no_more_missing_observations = 1;
end

% Compute the empirical mean of the observed variables.
dataset_info.descriptive.mean = nanmean(dataset_.data,1);

% Compute the empirical covariance matrix of the observed variables.
dataset_info.descriptive.covariance = nancovariance(dataset_.data);

% Compute the empirical correlation matrix of the observed variables.
normalization_matrix = diag(1./sqrt(diag(dataset_info.descriptive.covariance)));
dataset_info.descriptive.correlation = normalization_matrix*dataset_info.descriptive.covariance*normalization_matrix;

% Compute autocorrelation function.
dataset_info.descriptive.autocovariance = nanautocovariance(dataset_.data, options_.ar);

% Save raw data.
dataset_info.rawdata = dataset_.data;

% Prefilter the data if needed (remove the mean).
if isequal(options_.prefilter, 1)
    dataset_ = dataset_.detrend();
end
