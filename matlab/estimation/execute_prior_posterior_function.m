function oo_=execute_prior_posterior_function(posterior_function_name,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info,type)
%[oo_] = execute_prior_posterior_function(posterior_function_name,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info,type)
% This function executes a given function on draws of the posterior or prior distribution
%
% INPUTS
%   functionhandle               Handle to the function to be executed
%   M_           [structure]     Matlab/Octave structure describing the Model (initialized by dynare, see @ref{M_}).
%   options_     [structure]     Matlab/Octave structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_          [structure]     Matlab/Octave structure gathering the results (initialized by dynare, see @ref{oo_}).
%   estim_params_[structure]     Matlab/Octave structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   bayestopt_   [structure]     Matlab/Octave structure describing the parameter options (initialized by dynare, see @ref{bayestopt_}).
%   dataset_     [structure]     Matlab/Octave structure storing the dataset
%   dataset_info [structure]     Matlab/Octave structure storing the information about the dataset
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   oo_          [structure]     Matlab/Octave structure gathering the results (initialized by dynare, see @ref{oo_}).

% Copyright © 2013-2023 Dynare Team
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

[~,basename,extension] = fileparts(posterior_function_name);
if isempty(extension)
    extension = '.m';
end
fullname = [basename extension];
if ~strcmp(extension,'.m') %if not m-file
    error('The Posterior Function is not an m-file.')
elseif ~exist(fullname,'file') %if m-file, but does not exist
    error(['The Posterior Function ', fullname ,' was not found. Check the spelling.']);
end
%Create function handle
functionhandle=str2func(posterior_function_name);

n_draws=options_.sampling_draws;

if strcmpi(type,'posterior')
    % Get informations about the _posterior_draws files.
    % discard first mh_drop percent of the draws:
    CutSample(M_, options_, 'prior_posterior_function');
    % initialize metropolis draws
    options_.sub_draws = n_draws; % set draws for sampling; changed value is not returned to base workspace
    [error_flag, ~, options_] = metropolis_draw(1, options_, estim_params_, M_);
    if error_flag
        error('EXECUTE_POSTERIOR_FUNCTION: The draws could not be initialized')
    end
    n_draws = options_.sub_draws;
elseif strcmpi(type,'prior')
    % Get informations about the prior distribution.
    if isempty(bayestopt_)
        if ~isempty(estim_params_) && ~(isfield(estim_params_,'nvx') && (size(estim_params_.var_exo,1)+size(estim_params_.var_endo,1)+size(estim_params_.corrx,1)+size(estim_params_.corrn,1)+size(estim_params_.param_vals,1))==0)
            [~,estim_params_,bayestopt_,~,~,M_] = set_prior(estim_params_,M_,options_);
        else
            error('The prior distributions are not properly set up.')
        end
    end
    if exist([M_.fname '_prior_restrictions.m'])
        warning('prior_function currently does not support endogenous prior restrictions. They will be ignored. Consider using a posterior_function with nobs=1.')
    end
    Prior = dprior(bayestopt_, options_.prior_trunc);
else
    error('EXECUTE_POSTERIOR_FUNCTION: Unknown type!')
end

if strcmpi(type, 'prior')
    parameter_mat = Prior.draws(n_draws);
else
    parameter_mat = NaN(length(bayestopt_.p6), n_draws);
    for i = 1:n_draws
        parameter_mat(:,i) = GetOneDraw(type, M_, estim_params_, oo_, options_, bayestopt_);
    end
end

% Get output size
try
    junk = functionhandle(parameter_mat(:,1), M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
catch err
    fprintf('\nEXECUTE_POSTERIOR_FUNCTION: Execution of prior/posterior function led to an error. Execution cancelled.\n')
    rethrow(err)
end

% Initialize cell with number of columns
results_cell = cell(n_draws, columns(junk));

% Evaluate function on each draw
for i = 1:n_draws
    M_ = set_all_parameters(parameter_mat(:,i), estim_params_, M_);
    [results_cell(i,:)] = functionhandle(parameter_mat(:,i), M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
end

% Save results under oo_
oo_.(sprintf('%s_function_results', type)) = results_cell;
