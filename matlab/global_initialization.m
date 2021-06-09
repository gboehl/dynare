function global_initialization()
%function global_initialization()
% initializes global variables and options for DYNARE
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2018 Dynare Team
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

global oo_ M_ options_ estim_params_ bayestopt_ estimation_info ex0_ ys0_ dataset_ dataset_info
estim_params_ = [];
bayestopt_ = [];
dataset_=[];
dataset_info=[];

M_.dname = M_.fname;
M_.bvar = [];

estimation_info.empty_prior = struct(...
    'domain', [], 'interval', [], 'mean', [], ...
    'median', [], 'mode', [], 'shape', [], ...
    'shift', [], 'stdev', [], 'truncate', [], 'variance', []);
estimation_info.empty_options = struct(...
    'bounds',[], 'init', [], 'jscale', []);
estimation_info.subsamples.range = struct('date1', [], 'date2', []);
estimation_info.parameter.prior = estimation_info.empty_prior;
estimation_info.parameter.subsample_prior = estimation_info.empty_prior;
estimation_info.parameter.options = estimation_info.empty_options;
estimation_info.parameter.subsample_options = estimation_info.empty_options;
estimation_info.structural_innovation.prior = estimation_info.empty_prior;
estimation_info.structural_innovation.subsample_prior = estimation_info.empty_prior;
estimation_info.structural_innovation.options = estimation_info.empty_options;
estimation_info.structural_innovation.subsample_options = estimation_info.empty_options;
estimation_info.structural_innovation_corr.prior = estimation_info.empty_prior;
estimation_info.structural_innovation_corr.subsample_prior = estimation_info.empty_prior;
estimation_info.structural_innovation_corr.options = estimation_info.empty_options;
estimation_info.structural_innovation_corr.subsample_options = estimation_info.empty_options;
estimation_info.measurement_error.prior = estimation_info.empty_prior;
estimation_info.measurement_error.subsample_prior = estimation_info.empty_prior;
estimation_info.measurement_error.options = estimation_info.empty_options;
estimation_info.measurement_error.subsample_options = estimation_info.empty_options;
estimation_info.measurement_error_corr.prior = estimation_info.empty_prior;
estimation_info.measurement_error_corr.subsample_prior = estimation_info.empty_prior;
estimation_info.measurement_error_corr.options = estimation_info.empty_options;
estimation_info.measurement_error_corr.subsample_options = estimation_info.empty_options;
estimation_info.subsamples_index = {};
estimation_info.subsamples.range_index = {};
estimation_info.parameter_prior_index = {};
estimation_info.parameter_options_index = {};
estimation_info.parameter.range_index = {};
estimation_info.measurement_error_prior_index = {};
estimation_info.measurement_error_options_index = {};
estimation_info.measurement_error.range_index = {};
estimation_info.structural_innovation_prior_index = {};
estimation_info.structural_innovation_options_index = {};
estimation_info.structural_innovation.range_index = {};
estimation_info.measurement_error_corr_prior_index = {};
estimation_info.measurement_error_corr_options_index = {};
estimation_info.measurement_error_corr.range_index = {};
estimation_info.structural_innovation_corr_prior_index = {};
estimation_info.structural_innovation_corr_options_index = {};
estimation_info.structural_innovation_corr.range_index = {};
estimation_info.joint_parameter_prior_index = {};
estimation_info.joint_parameter = {'index','domain','interval','mean','median','mode','shape','shift','stdev','truncate','variance'};

oo_.exo_simul = [];
oo_.endo_simul = [];
ys0_ = [];
ex0_ = [];
oo_.dr = [];
oo_.exo_steady_state = [];
oo_.exo_det_steady_state = [];
oo_.exo_det_simul = [];
oo_.initval_series = dseries();

oo_.gui.ran_estimation = false;
oo_.gui.ran_stoch_simul = false;
oo_.gui.ran_calib_smoother = false;
oo_.gui.ran_perfect_foresight = false;
oo_.gui.ran_shock_decomposition = false;
oo_.gui.ran_realtime_shock_decomposition = false;

M_.params = [];
M_.endo_histval = [];
M_.exo_histval = [];
M_.exo_det_histval = [];
M_.Correlation_matrix = [];
M_.Correlation_matrix_ME = [];
M_.parameter_used_with_lead_lag = false;

M_.xref1.param = {};
M_.xref1.endo = {};
M_.xref1.exo = {};
M_.xref1.exo_det = {};

M_.xref2.param = {};
M_.xref2.endo = {};
M_.xref2.exo = {};
M_.xref2.exo_det = {};

M_.osr.param_names={};
M_.osr.param_indices=[];
M_.osr.param_bounds=[];
M_.osr.variable_weights=[];
M_.osr.variable_indices =[];

M_.instr_id=[];

% Set default options_ but keep global_init_file field if defined in the driver.
if isstruct(options_) && isfield(options_, 'global_init_file')
    global_init_file = options_.global_init_file;
    options_ = default_option_values(M_);
    options_.global_init_file = global_init_file;
else
    options_ = default_option_values(M_);
end

% initialize persistent variables in priordens()
priordens([],[],[],[],[],[],1);
% initialize persistent variables in dyn_first_order_solver()
dyn_first_order_solver();

% Set dynare random generator and seed.
set_dynare_seed('default');

% Load user configuration file.
if isfield(options_, 'global_init_file')
    if isfile(options_.global_init_file)
        try
            run(options_.global_init_file);
        catch
            error('Cannot evaluate global initialization file (%s)', options_.global_init_file)
        end
    else
        error('Cannot find global initialization file (%s).', options_.global_init_file)
    end
end

end

