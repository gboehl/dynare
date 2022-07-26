function [pdraws, STO_REDUCEDFORM, STO_MOMENTS, STO_DYNAMIC, STO_si_dDYNAMIC, STO_si_dREDUCEDFORM, STO_si_dMOMENTS, STO_dSPECTRUM, STO_dMINIMAL] = dynare_identification(options_ident, pdraws0)
%function [pdraws, STO_REDUCEDFORM, STO_MOMENTS, STO_DYNAMIC, STO_si_dDYNAMIC, STO_si_dREDUCEDFORM, STO_si_dMOMENTS, STO_dSPECTRUM, STO_dMINIMAL] = dynare_identification(options_ident, pdraws0)
% -------------------------------------------------------------------------
% This function is called, when the user specifies identification(...); in the mod file. It prepares all identification analysis:
% (1) set options, local and persistent variables for a new identification
%     analysis either for a single point or a MC Sample. It also displays and plots the results
% (2) load, display and plot a previously saved identification analysis
%
% Note 1: This function does not output the arguments to the workspace, but saves results to the folder identification
% Note 2: If you want to use this function directly in the mod file and workspace, you still have
%         to put identification in your mod file, otherwise the preprocessor won't provide all necessary objects
% =========================================================================
% INPUTS
%    * options_ident    [structure] identification options
%    * pdraws0          [SampleSize by totparam_nbr] optional: matrix of MC sample of model parameters
% -------------------------------------------------------------------------
% OUTPUTS
%    * pdraws               [matrix] MC sample of model params used
%    * STO_REDUCEDFORM,     [matrix] MC sample of entries in steady state and reduced form model solution (stacked vertically)
%    * STO_MOMENTS,         [matrix] MC sample of entries in theoretical first two moments (stacked vertically)
%    * STO_DYNAMIC,         [matrix] MC sample of entries in steady state and dynamic model derivatives (stacked vertically)
%    * STO_si_dDYNAMIC,     [matrix] MC sample of derivatives of steady state and dynamic derivatives
%    * STO_si_dREDUCEDFORM, [matrix] MC sample of derivatives of steady state and reduced form model solution
%    * STO_si_dMOMENTS      [matrix] MC sample of Iskrev (2010)'s J matrix
%    * STO_dSPECTRUM        [matrix] MC sample of Qu and Tkachenko (2012)'s \bar{G} matrix
%    * STO_dMINIMAL         [matrix] MC sample of Komunjer and Ng (2011)'s \bar{\Delta} matrix
% -------------------------------------------------------------------------
% This function is called by
%    * driver.m
%    * map_ident_.m
% -------------------------------------------------------------------------
% This function calls
%    * checkpath
%    * disp_identification
%    * dyn_waitbar
%    * dyn_waitbar_close
%    * get_all_parameters
%    * get_posterior_parameters
%    * get_the_name
%    * identification_analysis
%    * isoctave
%    * plot_identification
%    * prior_draw
%    * set_default_option
%    * set_prior
%    * skipline
%    * vnorm
% =========================================================================
% Copyright (C) 2010-2022 Dynare Team
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
% =========================================================================

global M_ options_ oo_ bayestopt_ estim_params_

% The TeX option crashes MATLAB R2014a run with "-nodisplay" option
% (as is done from the testsuite).
% Since we can’t directly test whether "-nodisplay" has been passed,
% we test for the "TOP_TEST_DIR" environment variable, which is set
% by the testsuite.
% Note that it was not tested whether the crash happens with more
% recent MATLAB versions, so when OLD_MATLAB_VERSION is increased,
% one should make a test before removing this workaround.
if options_.TeX && ~isoctave && matlab_ver_less_than('8.4') && ~isempty(getenv('TOP_TEST_DIR'))
    warning('Disabling TeX option due to a bug in MATLAB R2014a with -nodisplay')
    options_.TeX = false;
end

store_options_ = options_; % store options to restore them at the end
fname = M_.fname; %model name
dname = M_.dname; %model name

%turn warnings off, either globally or only relevant ids
if isoctave
    %warning('off'),
    warning('off','Octave:singular-matrix');
    warning('off','Octave:nearly-singular-matrix');
    warning('off','Octave:neg-dim-as-zero');
    warning('off','Octave:array-as-logical');
else
    %warning off;
    warning('off','MATLAB:rankDeficientMatrix');
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
    warning('off','MATLAB:specgraph:private:specgraph:UsingOnlyRealComponentOfComplexData');
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
    warning('off','MATLAB:divideByZero');
    warning('off','MATLAB:log:logOfZero');
end

%% Set all options in options_ident and create objects
options_ident = set_default_option(options_ident,'gsa_sample_file',0);
    % 0: do not use sample file
    % 1: triggers gsa prior sample
    % 2: triggers gsa Monte-Carlo sample (i.e. loads a sample corresponding to pprior=0 and ppost=0 in dynare_sensitivity options)
    % FILENAME: use sample file in provided path
options_ident = set_default_option(options_ident,'parameter_set','prior_mean');
    % 'calibration':      use values in M_.params and M_.Sigma_e to update estimated stderr, corr and model parameters (get_all_parameters)
    % 'prior_mode':       use values in bayestopt_.p5 prior to update estimated stderr, corr and model parameters
    % 'prior_mean':       use values in bayestopt_.p1 prior to update estimated stderr, corr and model parameters
    % 'posterior_mode':   use posterior mode values in estim_params_ to update estimated stderr, corr and model parameters (get_posterior_parameters('mode',...))
    % 'posterior_mean':   use posterior mean values in estim_params_ to update estimated stderr, corr and model parameters (get_posterior_parameters('mean',...))
    % 'posterior_median': use posterior median values in estim_params_ to update estimated stderr, corr and model parameters (get_posterior_parameters('median',...))
options_ident = set_default_option(options_ident,'load_ident_files',0);
    % 1: load previously computed analysis from identification/fname_identif.mat
options_ident = set_default_option(options_ident,'useautocorr',0);
    % 1: use autocorrelations in Iskrev (2010)'s theoretical second moments criteria
    % 0: use autocovariances in Iskrev (2010)'s theoretical second moments criteria
options_ident = set_default_option(options_ident,'ar',1);
    % number of lags to consider for autocovariances/autocorrelations in Iskrev (2010)'s criteria
options_ident = set_default_option(options_ident,'prior_mc',1);
    % size of Monte-Carlo sample of parameter draws
options_ident = set_default_option(options_ident,'prior_range',0);
    % 1: sample uniformly from prior ranges implied by the prior specifications (overwrites prior shape when prior_mc > 1)
    % 0: sample from specified prior distributions (when prior_mc > 1)
options_ident = set_default_option(options_ident,'periods',300);
    % length of stochastic simulation to compute simulated moment uncertainty (when analytic Hessian is not available)
options_ident = set_default_option(options_ident,'replic',100);
    % number of replicas to compute simulated moment uncertainty (when analytic Hessian is not available)
options_ident = set_default_option(options_ident,'advanced',0);
    % 1: show a more detailed analysis based on reduced-form model solution and dynamic model derivatives. Further, performs a brute force
    %    search of the groups of parameters best reproducing the behavior of each single parameter.
options_ident = set_default_option(options_ident,'normalize_jacobians',1);
    % 1: normalize Jacobians by either rescaling each row by its largest element in absolute value or for Gram (or Hessian-type) matrices by transforming into correlation-type matrices
options_ident = set_default_option(options_ident,'grid_nbr',5000);
    % number of grid points in [-pi;pi] to approximate the integral to compute Qu and Tkachenko (2012)'s criteria
    % note that grid_nbr needs to be even and actually we use (grid_nbr+1) points, as we add the 0 frequency and use symmetry
    if mod(options_ident.grid_nbr,2) ~= 0
        options_ident.grid_nbr = options_ident.grid_nbr+1;
        fprintf('IDENTIFICATION: ''grid_nbr'' needs to be even, hence it is reset to %d\n',options_ident.grid_nbr)
        if mod(options_ident.grid_nbr,2) ~= 0
            error('IDENTIFICATION: You need to set an even value for ''grid_nbr''');
        end
    end
options_ident = set_default_option(options_ident,'tol_rank','robust');
    % tolerance level used for rank computations
options_ident = set_default_option(options_ident,'tol_deriv',1.e-8);
    % tolerance level for selecting columns of non-zero derivatives
options_ident = set_default_option(options_ident,'tol_sv',1.e-3);
    % tolerance level for selecting non-zero singular values in identification_checks.m
options_ident = set_default_option(options_ident,'schur_vec_tol',1e-11);
    % tolerance level used to find nonstationary variables in Schur decomposition of the transition matrix.

%check whether to compute identification strength based on information matrix
if ~isfield(options_ident,'no_identification_strength')
    options_ident.no_identification_strength = 0;
else
    options_ident.no_identification_strength = 1;
end
%check whether to compute and display identification criteria based on steady state and reduced-form model solution
if ~isfield(options_ident,'no_identification_reducedform')
    options_ident.no_identification_reducedform = 0;
else
    options_ident.no_identification_reducedform = 1;
end
%check whether to compute and display identification criteria based on Iskrev (2010), i.e. derivative of first two moments
if ~isfield(options_ident,'no_identification_moments')
    options_ident.no_identification_moments = 0;
else
    options_ident.no_identification_moments = 1;
end
%check whether to compute and display identification criteria based on Komunjer and Ng (2011), i.e. derivative of first moment, minimal state space system and observational equivalent spectral density transformation
if ~isfield(options_ident,'no_identification_minimal')
    options_ident.no_identification_minimal = 0;
else
    options_ident.no_identification_minimal = 1;
end
%Check whether to compute and display identification criteria based on Qu and Tkachenko (2012), i.e. Gram matrix of derivatives of first moment plus outer product of derivatives of spectral density
if ~isfield(options_ident,'no_identification_spectrum')
    options_ident.no_identification_spectrum = 0;
else
    options_ident.no_identification_spectrum = 1;
end

%overwrite setting, as identification strength and advanced need both reducedform and moments criteria
if (isfield(options_ident,'no_identification_strength') &&  options_ident.no_identification_strength == 0) || (options_ident.advanced == 1)
    options_ident.no_identification_reducedform = 0;
    options_ident.no_identification_moments = 0;
end

%overwrite setting, as dynare_sensitivity does not make use of spectrum and minimal system
if isfield(options_,'opt_gsa') && isfield(options_.opt_gsa,'identification') && options_.opt_gsa.identification == 1
    options_ident.no_identification_minimal = 1;
    options_ident.no_identification_spectrum = 1;
end

%Deal with non-stationary cases
if isfield(options_ident,'diffuse_filter')
    options_ident.lik_init=3; %overwrites any other lik_init set
    options_.diffuse_filter=options_ident.diffuse_filter; %make options_ inherit diffuse filter; will overwrite any conflicting lik_init in dynare_estimation_init
else
    if options_.diffuse_filter==1 %warning if estimation with diffuse filter was done, but not passed
        fprintf('WARNING IDENTIFICATION: Previously the diffuse_filter option was used, but it was not passed to the identification command. This may result in problems if your model contains unit roots.\n');
    end
end
options_ident = set_default_option(options_ident,'lik_init',1);
options_.lik_init=options_ident.lik_init; %make options_ inherit lik_init
if options_ident.lik_init==3 %user specified diffuse filter using the lik_init option
    options_ident.analytic_derivation=0; %diffuse filter not compatible with analytic derivation
end
    % Type of initialization of Kalman filter:
    % 1: stationary models: initial matrix of variance of error of forecast is set equal to the unconditional variance of the state variables
    % 2: nonstationary models: wide prior is used with an initial matrix of variance of the error of forecast diagonal with 10 on the diagonal (follows the suggestion of Harvey and Phillips(1979))
    % 3: nonstationary models: use a diffuse filter (use rather the diffuse_filter option)
    % 4: filter is initialized with the fixed point of the Riccati equation
    % 5:  i) option 2 for non-stationary elements by setting their initial variance in the forecast error matrix to 10 on the diagonal and all co-variances to 0 and
    %    ii) option 1 for the stationary elements
options_ident = set_default_option(options_ident,'analytic_derivation',1);
    % 1: analytic derivation of gradient and hessian of likelihood in dsge_likelihood.m, only works for stationary models, i.e. kalman_algo<3
options_ident = set_default_option(options_ident,'order',1);
    % 1: first-order perturbation approximation, identification is based on linear state space system
    % 2: second-order perturbation approximation, identification is based on second-order pruned state space system
    % 3: third-order perturbation approximation, identification is based on third-order pruned state space system

%overwrite values in options_, note that options_ is restored at the end of the function
if isfield(options_ident,'prior_trunc')
    options_.prior_trunc=options_ident.prior_trunc;
        % sets truncation of prior
end
if isfield(options_ident,'TeX')
    options_.TeX=options_ident.TeX;
        % requests printing of results and graphs in LaTeX
end
if isfield(options_ident,'nograph')
    options_.nograph=options_ident.nograph;
        % do not display and do not save graphs
end
if isfield(options_ident,'nodisplay')
    options_.nodisplay=options_ident.nodisplay;
        % do not display, but save graphs
end
if isfield(options_ident,'graph_format')
    options_.graph_format=options_ident.graph_format;
        % specify file formats to save graphs: eps, pdf, fig, none (do not save but only display)
end

% check for external draws, i.e. set pdraws0 for a gsa analysis
if options_ident.gsa_sample_file
    GSAFolder = checkpath('gsa',dname);
    if options_ident.gsa_sample_file==1
        load([GSAFolder,filesep,fname,'_prior'],'lpmat','lpmat0','istable');
    elseif options_ident.gsa_sample_file==2
        load([GSAFolder,filesep,fname,'_mc'],'lpmat','lpmat0','istable');
    else
        load([GSAFolder,filesep,options_ident.gsa_sample_file],'lpmat','lpmat0','istable');
    end
    if isempty(istable)
        istable=1:size(lpmat,1);
    end
    if ~isempty(lpmat0)
        lpmatx=lpmat0(istable,:);
    else
        lpmatx=[];
    end
    pdraws0 = [lpmatx lpmat(istable,:)];
    clear lpmat lpmat0 istable;
elseif nargin==1
    pdraws0=[];
end
external_sample=0;
if nargin==2 || ~isempty(pdraws0)
    % change settings if there is an external sample provided as input argument
    options_ident.prior_mc = size(pdraws0,1);
    options_ident.load_ident_files = 0;
    external_sample = 1;
end

% Check if estimated_params block is provided
if isempty(estim_params_)
    prior_exist = 0;
    %reset some options
    options_ident.prior_mc = 1;
    options_ident.prior_range = 0;
    options_.identification_check_endogenous_params_with_no_prior = true; %needed to trigger endogenous steady state parameter check in dynare_estimation_init
else
    prior_exist = 1;
    parameters = options_ident.parameter_set;
end

% overwrite settings in options_ and prepare to call dynare_estimation_init
options_.order = options_ident.order;
if options_ident.order > 1
    %order>1 is not compatible with analytic derivation in dsge_likelihood.m
    options_ident.analytic_derivation=0;
    %order>1 is based on pruned state space system
    options_.pruning = true;
end
if options_ident.order == 3
    options_.k_order_solver = 1;
end
options_.ar = options_ident.ar;
options_.prior_mc = options_ident.prior_mc;
options_.schur_vec_tol = options_ident.schur_vec_tol;
options_.nomoments = 0;
options_.analytic_derivation=options_ident.analytic_derivation;
    % 1: analytic derivation of gradient and hessian of likelihood in dsge_likelihood.m, only works for stationary models, i.e. kalman_algo<3
options_ = set_default_option(options_,'datafile','');
options_.mode_compute = 0;
options_.plot_priors = 0;
options_.smoother = 1;
options_.options_ident = [];
[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = dynare_estimation_init(M_.endo_names, fname, 1, M_, options_, oo_, estim_params_, bayestopt_);

% set method to compute identification Jacobians (kronflag). Default:0
options_ident = set_default_option(options_ident,'analytic_derivation_mode', options_.analytic_derivation_mode); % if not set by user, inherit default global one
    %  0: efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2012)
    %  1: kronecker products method to compute analytical derivatives as in Iskrev (2010) (only for order=1)
    % -1: numerical two-sided finite difference method to compute numerical derivatives of all identification Jacobians using function identification_numerical_objective.m (previously thet2tau.m)
    % -2: numerical two-sided finite difference method to compute numerically dYss, dg1, dg2, dg3, d2Yss and d2g1, the identification Jacobians are then computed analytically as with 0

if options_.discretionary_policy || options_.ramsey_policy
    if options_ident.analytic_derivation_mode~=-1
        fprintf('dynare_identification: discretionary_policy and ramsey_policy require analytic_derivation_mode=-1. Resetting the option.')
        options_ident.analytic_derivation_mode=-1;
    end
end
    
options_.analytic_derivation_mode = options_ident.analytic_derivation_mode; %overwrite setting in options_

% initialize persistent variables in prior_draw
if prior_exist
    if any(bayestopt_.pshape > 0)
        if options_ident.prior_range
            %sample uniformly from prior ranges (overwrite prior specification)
            prior_draw(bayestopt_, options_.prior_trunc, true);
        else
            %sample from prior distributions
            prior_draw(bayestopt_, options_.prior_trunc, false);
        end
    else
        options_ident.prior_mc = 1; %only one single point
    end
end

% check if identification directory is already created
IdentifDirectoryName = CheckPath('identification',dname);

% Create indices for stderr, corr and model parameters
if prior_exist % use estimated_params block
    indpmodel = []; %initialize index for model parameters
    if ~isempty(estim_params_.param_vals)
        indpmodel = estim_params_.param_vals(:,1); %values correspond to parameters declaration order, row number corresponds to order in estimated_params
    end
    indpstderr=[]; %initialize index for stderr parameters
    if ~isempty(estim_params_.var_exo)
        indpstderr = estim_params_.var_exo(:,1); %values correspond to varexo declaration order, row number corresponds to order in estimated_params
    end
    indpcorr=[]; %initialize matrix for corr paramters
    if ~isempty(estim_params_.corrx)
        indpcorr = estim_params_.corrx(:,1:2); %values correspond to varexo declaration order, row number corresponds to order in estimated_params
    end
    totparam_nbr = length(bayestopt_.name); % number of estimated stderr, corr and model parameters as declared in estimated_params
    modparam_nbr = estim_params_.np;        % number of model parameters as declared in estimated_params
    stderrparam_nbr = estim_params_.nvx;    % number of stderr parameters
    corrparam_nbr = estim_params_.ncx;      % number of corr parameters
    if estim_params_.nvn || estim_params_.ncn %nvn is number of stderr parameters and ncn is number of corr parameters of measurement innovations as declared in estimated_params
        error('Identification does not (yet) support measurement errors. Instead, define them explicitly as varexo and provide measurement equations in the model definition.')
    end
    name = cell(totparam_nbr,1);     %initialize cell for parameter names
    name_tex = cell(totparam_nbr,1); %initialize cell for TeX parameter names
    for jj=1:totparam_nbr
        if options_.TeX
            [param_name_temp, param_name_tex_temp]= get_the_name(jj,options_.TeX,M_,estim_params_,options_);
            name_tex{jj,1} = strrep(param_name_tex_temp,'$',''); %ordering corresponds to estimated_params
            name{jj,1} = param_name_temp; %ordering corresponds to estimated_params
        else
            param_name_temp = get_the_name(jj,options_.TeX,M_,estim_params_,options_);
            name{jj,1} = param_name_temp; %ordering corresponds to estimated_params
        end
    end
else % no estimated_params block, choose all model parameters and all stderr parameters, but no corr parameters
    indpmodel = 1:M_.param_nbr; %all model parameters
    indpstderr = 1:M_.exo_nbr;  %all stderr parameters
    indpcorr = [];              %no corr parameters
    stderrparam_nbr = M_.exo_nbr;
    corrparam_nbr = 0;
    modparam_nbr = M_.param_nbr;
    totparam_nbr = modparam_nbr+stderrparam_nbr;
    name = cellfun(@(x) horzcat('SE_', x), M_.exo_names, 'UniformOutput', false); %names for stderr parameters
    name = vertcat(name, M_.param_names);
    name_tex = cellfun(@(x) horzcat('$ SE_{', x, '} $'), M_.exo_names, 'UniformOutput', false);
    name_tex = vertcat(name_tex, M_.param_names_tex);
    if ~isequal(M_.H,0)
        fprintf('\ndynare_identification:: Identification does not support measurement errors (yet) and will ignore them in the following. To test their identifiability, instead define them explicitly as varexo and provide measurement equations in the model definition.\n')
    end
end
options_ident.name_tex = name_tex;

fprintf('\n======== Identification Analysis ========\n')
if options_ident.order > 1
    fprintf('Based on Pruned State Space System (order=%d)\n',options_ident.order);
end
skipline()
if totparam_nbr < 2
    options_ident.advanced = 0;
    fprintf('There is only one parameter to study for identitification. The advanced option is re-set to 0.\n')
end

% settings dependent on number of parameters
options_ident = set_default_option(options_ident,'max_dim_cova_group',min([2,totparam_nbr-1]));
options_ident.max_dim_cova_group = min([options_ident.max_dim_cova_group,totparam_nbr-1]);
    % In brute force search (ident_bruteforce.m) when advanced=1 this option sets the maximum dimension of groups of parameters that best reproduce the behavior of each single model parameter

options_ident = set_default_option(options_ident,'checks_via_subsets',0);
    % 1: uses identification_checks_via_subsets.m to compute problematic parameter combinations
    % 0: uses identification_checks.m to compute problematic parameter combinations [default]
options_ident = set_default_option(options_ident,'max_dim_subsets_groups',min([4,totparam_nbr-1]));
    % In identification_checks_via_subsets.m, when checks_via_subsets=1, this option sets the maximum dimension of groups of parameters for which the corresponding rank criteria is checked


% store identification options
options_.options_ident = options_ident;
store_options_ident = options_ident;

% get some options for quick reference
iload = options_ident.load_ident_files;
SampleSize = options_ident.prior_mc;

if iload <=0
    %% Perform new identification analysis, i.e. do not load previous analysis
    if prior_exist
        % use information from estimated_params block
        params = set_prior(estim_params_,M_,options_)';
        if all(bayestopt_.pshape == 0)
            % only bounds are specified in estimated_params
            parameters = 'ML_Starting_value';
            parameters_TeX = 'ML starting value';
            fprintf('Testing ML Starting value\n');
        else
            % use user-defined option
            switch parameters
                case 'calibration'
                    parameters_TeX = 'Calibration';
                    fprintf('Testing calibration\n');
                    params(1,:) = get_all_parameters(estim_params_,M_);
                case 'posterior_mode'
                    parameters_TeX = 'Posterior mode';
                    fprintf('Testing posterior mode\n');
                    params(1,:) = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
                case 'posterior_mean'
                    parameters_TeX = 'Posterior mean';
                    fprintf('Testing posterior mean\n');
                    params(1,:) = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
                case 'posterior_median'
                    parameters_TeX = 'Posterior median';
                    fprintf('Testing posterior median\n');
                    params(1,:) = get_posterior_parameters('median',M_,estim_params_,oo_,options_);
                case 'prior_mode'
                    parameters_TeX = 'Prior mode';
                    fprintf('Testing prior mode\n');
                    params(1,:) = bayestopt_.p5(:);
                case 'prior_mean'
                    parameters_TeX = 'Prior mean';
                    fprintf('Testing prior mean\n');
                    params(1,:) = bayestopt_.p1;
                otherwise
                    fprintf('The option parameter_set has to be equal to: ''calibration'', ''posterior_mode'', ''posterior_mean'', ''posterior_median'', ''prior_mode'' or ''prior_mean''.\n');
                    error('IDENTIFICATION: The option ''parameter_set'' has an invalid value');
            end
        end
    else
        % no estimated_params block is available, all stderr and model parameters, but no corr parameters are chosen
        params = [sqrt(diag(M_.Sigma_e))', M_.params'];
        parameters = 'Current_params';
        parameters_TeX = 'Current parameter values';
        fprintf('Testing all current stderr and model parameter values\n');
    end
    options_ident.tittxt = parameters; %title text for graphs and figures
    % perform identification analysis for single point
    [ide_moments_point, ide_spectrum_point, ide_minimal_point, ide_hess_point, ide_reducedform_point, ide_dynamic_point, derivatives_info_point, info, error_indicator_point] = ...
        identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end implies initialization of persistent variables
    if info(1)~=0
        % there are errors in the solution algorithm
        message = get_error_message(info,options_);
        fprintf('-----------\n');
        fprintf('The model does not solve for %s (info = %d: %s)\n', parameters, info(1), message);
        fprintf('-----------\n');
        if any(bayestopt_.pshape)
            % if there are errors in the solution algorithm, try to sample a different point from the prior
            fprintf('Try sampling up to 50 parameter sets from the prior.\n');
            kk=0;
            while kk<50 && info(1)
                kk=kk+1;
                params = prior_draw();
                options_ident.tittxt = 'Random_prior_params'; %title text for graphs and figures
                % perform identification analysis
                [ide_moments_point, ide_spectrum_point, ide_minimal_point, ide_hess_point, ide_reducedform_point, ide_dynamic_point, derivatives_info_point, info, error_indicator_point] = ...
                    identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1);
            end
        end
        if info(1)
            fprintf('\n-----------\n');
            fprintf('Identification stopped:\n');
            if any(bayestopt_.pshape)
                fprintf('The model did not solve for any of 50 attempts of random samples from the prior\n');
            end
            fprintf('-----------\n');
            return
        else
            % found a (random) point that solves the model
            fprintf('Found a random draw from the priors that solves the model:\n');
            disp(params);
            fprintf('Identification now continues for this draw.');
            parameters = 'Random_prior_params';
            parameters_TeX = 'Random prior parameter draw';
        end
    end
    ide_hess_point.params = params;
    % save all output into identification folder
    save([IdentifDirectoryName '/' fname                '_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point', 'ide_hess_point', 'ide_reducedform_point', 'ide_dynamic_point', 'store_options_ident');
    save([IdentifDirectoryName '/' fname '_' parameters '_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point', 'ide_hess_point', 'ide_reducedform_point', 'ide_dynamic_point', 'store_options_ident');
    % display results of identification analysis
    disp_identification(params, ide_reducedform_point, ide_moments_point, ide_spectrum_point, ide_minimal_point, name, options_ident);
    if ~options_ident.no_identification_strength && ~options_.nograph && ~error_indicator_point.identification_strength && ~error_indicator_point.identification_moments
        % plot (i) identification strength and sensitivity measure based on the moment information matrix and (ii) plot advanced analysis graphs
        plot_identification(params, ide_moments_point, ide_hess_point, ide_reducedform_point, ide_dynamic_point, options_ident.advanced, parameters, name, IdentifDirectoryName, parameters_TeX, name_tex);
    end

    if SampleSize > 1
        % initializations for Monte Carlo Analysis
        fprintf('\nMonte Carlo Testing\n');
        h = dyn_waitbar(0,'Monte Carlo identification checks ...');
        iteration  = 0; % initialize counter for admissable draws
        run_index  = 0; % initialize counter for admissable draws after saving previous draws to file(s)
        file_index = 0; % initialize counter for files (if options_.MaxNumberOfBytes is reached, we store results in files)
        options_MC = options_ident; %store options structure for Monte Carlo analysis
        options_MC.advanced = 0;    %do not run advanced checking in a Monte Carlo analysis
        options_ident.checks_via_subsets = 0; % for Monte Carlo analysis currently only identification_checks and not identification_checks_via_subsets is supported
    else
        iteration = 1; % iteration equals SampleSize and we are finished
        pdraws = [];   % to have output object otherwise map_ident.m may crash
    end
    while iteration < SampleSize
        if external_sample
            params = pdraws0(iteration+1,:); % loaded draws
        else
            params = prior_draw(); % new random draw from prior
        end
        options_ident.tittxt = []; % clear title text for graphs and figures
        % run identification analysis
        [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_dynamic, ide_derivatives_info, info, error_indicator] = ...
            identification_analysis(params, indpmodel, indpstderr, indpcorr, options_MC, dataset_info, prior_exist, 0); % the 0 implies that we do not initialize persistent variables anymore

        if iteration==0 && info(1)==0 % preallocate storage in the first admissable run
            delete([IdentifDirectoryName '/' fname '_identif_*.mat']) % delete previously saved results
            MAX_RUNS_BEFORE_SAVE_TO_FILE = min(SampleSize,ceil(options_.MaxNumberOfBytes/(size(ide_reducedform.si_dREDUCEDFORM,1)*totparam_nbr)/8)); % set how many runs can be stored before we save to files
            pdraws = zeros(SampleSize,totparam_nbr); % preallocate storage for draws in each row

            % preallocate storage for dynamic model
            STO_si_dDYNAMIC          = zeros([size(ide_dynamic.si_dDYNAMIC, 1), modparam_nbr, MAX_RUNS_BEFORE_SAVE_TO_FILE]);
            STO_DYNAMIC              = zeros(size(ide_dynamic.DYNAMIC, 1), SampleSize);
            IDE_DYNAMIC.ind_dDYNAMIC = ide_dynamic.ind_dDYNAMIC;
            IDE_DYNAMIC.ino          = zeros(SampleSize, modparam_nbr);
            IDE_DYNAMIC.ind0         = zeros(SampleSize, modparam_nbr);
            IDE_DYNAMIC.jweak        = zeros(SampleSize, modparam_nbr);
            IDE_DYNAMIC.jweak_pair   = zeros(SampleSize, modparam_nbr*(modparam_nbr+1)/2);
            IDE_DYNAMIC.cond         = zeros(SampleSize, 1);
            IDE_DYNAMIC.Mco          = zeros(SampleSize, modparam_nbr);

            % preallocate storage for reduced form
            if ~options_MC.no_identification_reducedform
                STO_si_dREDUCEDFORM              = zeros([size(ide_reducedform.si_dREDUCEDFORM, 1), totparam_nbr, MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                STO_REDUCEDFORM                  = zeros(size(ide_reducedform.REDUCEDFORM, 1), SampleSize);
                IDE_REDUCEDFORM.ind_dREDUCEDFORM = ide_reducedform.ind_dREDUCEDFORM;
                IDE_REDUCEDFORM.ino              = zeros(SampleSize, 1);
                IDE_REDUCEDFORM.ind0             = zeros(SampleSize, totparam_nbr);
                IDE_REDUCEDFORM.jweak            = zeros(SampleSize, totparam_nbr);
                IDE_REDUCEDFORM.jweak_pair       = zeros(SampleSize, totparam_nbr*(totparam_nbr+1)/2);
                IDE_REDUCEDFORM.cond             = zeros(SampleSize, 1);
                IDE_REDUCEDFORM.Mco              = zeros(SampleSize, totparam_nbr);
            else
                STO_si_dREDUCEDFORM = {};
                STO_REDUCEDFORM     = {};
                IDE_REDUCEDFORM     = {};
            end

            % preallocate storage for moments
            if ~options_MC.no_identification_moments && ~isempty(fieldnames(ide_moments))
                STO_si_dMOMENTS          = zeros([size(ide_moments.si_dMOMENTS, 1), totparam_nbr, MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                STO_MOMENTS              = zeros(size(ide_moments.MOMENTS, 1), SampleSize);
                IDE_MOMENTS.ind_dMOMENTS = ide_moments.ind_dMOMENTS;
                IDE_MOMENTS.ino          = zeros(SampleSize, 1);
                IDE_MOMENTS.ind0         = zeros(SampleSize, totparam_nbr);
                IDE_MOMENTS.jweak        = zeros(SampleSize, totparam_nbr);
                IDE_MOMENTS.jweak_pair   = zeros(SampleSize, totparam_nbr*(totparam_nbr+1)/2);
                IDE_MOMENTS.cond         = zeros(SampleSize, 1);
                IDE_MOMENTS.Mco          = zeros(SampleSize, totparam_nbr);
                IDE_MOMENTS.S            = zeros(SampleSize, min(8, totparam_nbr));
                IDE_MOMENTS.V            = zeros(SampleSize, totparam_nbr, min(8, totparam_nbr));
            else
                STO_si_dMOMENTS = {};
                STO_MOMENTS     = {};
                IDE_MOMENTS     = {};
            end

            % preallocate storage for spectrum
            if ~options_MC.no_identification_spectrum && ~isempty(fieldnames(ide_spectrum))
                STO_dSPECTRUM              = zeros([size(ide_spectrum.dSPECTRUM, 1), size(ide_spectrum.dSPECTRUM, 2), MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                IDE_SPECTRUM.ind_dSPECTRUM = ide_spectrum.ind_dSPECTRUM;
                IDE_SPECTRUM.ino           = zeros(SampleSize, 1);
                IDE_SPECTRUM.ind0          = zeros(SampleSize, totparam_nbr);
                IDE_SPECTRUM.jweak         = zeros(SampleSize, totparam_nbr);
                IDE_SPECTRUM.jweak_pair    = zeros(SampleSize, totparam_nbr*(totparam_nbr+1)/2);
                IDE_SPECTRUM.cond          = zeros(SampleSize, 1);
                IDE_SPECTRUM.Mco           = zeros(SampleSize, totparam_nbr);
            else
                STO_dSPECTRUM = {};
                IDE_SPECTRUM  = {};
            end

            % preallocate storage for minimal system
            if ~options_MC.no_identification_minimal && ~isempty(fieldnames(ide_minimal)) && ide_minimal.minimal_state_space==1
                STO_dMINIMAL             = zeros([size(ide_minimal.dMINIMAL, 1), size(ide_minimal.dMINIMAL, 2), MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                IDE_MINIMAL.ind_dMINIMAL = ide_minimal.ind_dMINIMAL;
                IDE_MINIMAL.ino          = zeros(SampleSize, 1);
                IDE_MINIMAL.ind0         = zeros(SampleSize, totparam_nbr);
                IDE_MINIMAL.jweak        = zeros(SampleSize, totparam_nbr);
                IDE_MINIMAL.jweak_pair   = zeros(SampleSize, totparam_nbr*(totparam_nbr+1)/2);
                IDE_MINIMAL.cond         = zeros(SampleSize, 1);
                IDE_MINIMAL.Mco          = zeros(SampleSize, totparam_nbr);
                IDE_MINIMAL.minimal_state_space = zeros(SampleSize, 1);
            else
                STO_dMINIMAL = {};
                IDE_MINIMAL  = {};
            end
        end

        if info(1)==0 % if admissable draw
            iteration = iteration + 1; %increase total index of admissable draws
            run_index = run_index + 1; %increase index of admissable draws after saving to files
            pdraws(iteration,:) = params; % store draw

            % store results for steady state and dynamic model derivatives
            STO_DYNAMIC(:,iteration)            = ide_dynamic.DYNAMIC;
            STO_si_dDYNAMIC(:,:,run_index)      = ide_dynamic.si_dDYNAMIC;
            IDE_DYNAMIC.cond(iteration,1)       = ide_dynamic.cond;
            IDE_DYNAMIC.ino(iteration,1)        = ide_dynamic.ino;
            IDE_DYNAMIC.ind0(iteration,:)       = ide_dynamic.ind0;
            IDE_DYNAMIC.jweak(iteration,:)      = ide_dynamic.jweak;
            IDE_DYNAMIC.jweak_pair(iteration,:) = ide_dynamic.jweak_pair;
            IDE_DYNAMIC.Mco(iteration,:)        = ide_dynamic.Mco;

            % store results for reduced form
            if ~options_MC.no_identification_reducedform && ~error_indicator.identification_reducedform
                STO_REDUCEDFORM(:,iteration)            = ide_reducedform.REDUCEDFORM;
                STO_si_dREDUCEDFORM(:,:,run_index)      = ide_reducedform.si_dREDUCEDFORM;
                IDE_REDUCEDFORM.cond(iteration,1)       = ide_reducedform.cond;
                IDE_REDUCEDFORM.ino(iteration,1)        = ide_reducedform.ino;
                IDE_REDUCEDFORM.ind0(iteration,:)       = ide_reducedform.ind0;
                IDE_REDUCEDFORM.jweak(iteration,:)      = ide_reducedform.jweak;
                IDE_REDUCEDFORM.jweak_pair(iteration,:) = ide_reducedform.jweak_pair;
                IDE_REDUCEDFORM.Mco(iteration,:)        = ide_reducedform.Mco;
            end

            % store results for moments
            if ~options_MC.no_identification_moments && ~error_indicator.identification_moments
                STO_MOMENTS(:,iteration)            = ide_moments.MOMENTS;
                STO_si_dMOMENTS(:,:,run_index)      = ide_moments.si_dMOMENTS;
                IDE_MOMENTS.cond(iteration,1)       = ide_moments.cond;
                IDE_MOMENTS.ino(iteration,1)        = ide_moments.ino;
                IDE_MOMENTS.ind0(iteration,:)       = ide_moments.ind0;
                IDE_MOMENTS.jweak(iteration,:)      = ide_moments.jweak;
                IDE_MOMENTS.jweak_pair(iteration,:) = ide_moments.jweak_pair;
                IDE_MOMENTS.Mco(iteration,:)        = ide_moments.Mco;
                IDE_MOMENTS.S(iteration,:)          = ide_moments.S;
                IDE_MOMENTS.V(iteration,:,:)        = ide_moments.V;
            end

            % store results for spectrum
            if ~options_MC.no_identification_spectrum && ~error_indicator.identification_spectrum
                STO_dSPECTRUM(:,:,run_index)         = ide_spectrum.dSPECTRUM;
                IDE_SPECTRUM.cond(iteration,1)       = ide_spectrum.cond;
                IDE_SPECTRUM.ino(iteration,1)        = ide_spectrum.ino;
                IDE_SPECTRUM.ind0(iteration,:)       = ide_spectrum.ind0;
                IDE_SPECTRUM.jweak(iteration,:)      = ide_spectrum.jweak;
                IDE_SPECTRUM.jweak_pair(iteration,:) = ide_spectrum.jweak_pair;
                IDE_SPECTRUM.Mco(iteration,:)        = ide_spectrum.Mco;
            end

            % store results for minimal system
            if ~options_MC.no_identification_minimal 
                if ~error_indicator.identification_minimal
                    STO_dMINIMAL(:,:,run_index)                  = ide_minimal.dMINIMAL;
                    IDE_MINIMAL.cond(iteration,1)                = ide_minimal.cond;
                    IDE_MINIMAL.ino(iteration,1)                 = ide_minimal.ino;
                    IDE_MINIMAL.ind0(iteration,:)                = ide_minimal.ind0;
                    IDE_MINIMAL.jweak(iteration,:)               = ide_minimal.jweak;
                    IDE_MINIMAL.jweak_pair(iteration,:)          = ide_minimal.jweak_pair;
                    IDE_MINIMAL.Mco(iteration,:)                 = ide_minimal.Mco;
                    IDE_MINIMAL.minimal_state_space(iteration,:) = ide_minimal.minimal_state_space;
                else
                    IDE_MINIMAL.minimal_state_space(iteration,:) = ide_minimal.minimal_state_space;
                end
            end

            % save results to file: either to avoid running into memory issues, i.e. (run_index==MAX_RUNS_BEFORE_SAVE_TO_FILE) or if finished (iteration==SampleSize)
            if run_index==MAX_RUNS_BEFORE_SAVE_TO_FILE || iteration==SampleSize
                file_index = file_index + 1;
                if run_index < MAX_RUNS_BEFORE_SAVE_TO_FILE
                    %we are finished (iteration == SampleSize), so get rid of additional storage
                    STO_si_dDYNAMIC = STO_si_dDYNAMIC(:,:,1:run_index);
                    if ~options_MC.no_identification_reducedform
                        STO_si_dREDUCEDFORM = STO_si_dREDUCEDFORM(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_moments
                        STO_si_dMOMENTS = STO_si_dMOMENTS(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_spectrum
                        STO_dSPECTRUM = STO_dSPECTRUM(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_minimal
                        STO_dMINIMAL = STO_dMINIMAL(:,:,1:run_index);
                    end
                end
                save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_dDYNAMIC');
                STO_si_dDYNAMIC = zeros(size(STO_si_dDYNAMIC)); % reset storage
                if ~options_MC.no_identification_reducedform
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_dREDUCEDFORM', '-append');
                    STO_si_dREDUCEDFORM = zeros(size(STO_si_dREDUCEDFORM)); % reset storage
                end
                if ~options_MC.no_identification_moments
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_dMOMENTS','-append');
                    STO_si_dMOMENTS = zeros(size(STO_si_dMOMENTS)); % reset storage
                end
                if ~options_MC.no_identification_spectrum
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_dSPECTRUM','-append');
                    STO_dSPECTRUM = zeros(size(STO_dSPECTRUM)); % reset storage
                end
                if ~options_MC.no_identification_minimal
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_dMINIMAL','-append');
                    STO_dMINIMAL = zeros(size(STO_dMINIMAL)); % reset storage
                end
                run_index = 0; % reset index
            end
            if SampleSize > 1
                dyn_waitbar(iteration/SampleSize, h, ['MC identification checks ', int2str(iteration), '/', int2str(SampleSize)]);
            end
        end
    end

    if SampleSize > 1
        dyn_waitbar_close(h);
        normalize_STO_DYNAMIC = std(STO_DYNAMIC,0,2);
        if ~options_MC.no_identification_reducedform
            normalize_STO_REDUCEDFORM = std(STO_REDUCEDFORM,0,2);
        end
        if ~options_MC.no_identification_moments
            normalize_STO_MOMENTS = std(STO_MOMENTS,0,2);
        end
        if ~options_MC.no_identification_minimal
            normalize_STO_MINIMAL = 1; %not used (yet)
        end
        if ~options_MC.no_identification_spectrum
            normalize_STO_SPECTRUM = 1; %not used (yet)
        end
        normaliz1 = std(pdraws);
        iter = 0;
        for ifile_index = 1:file_index
            load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_dDYNAMIC');
            maxrun_dDYNAMIC = size(STO_si_dDYNAMIC,3);
            if ~options_MC.no_identification_reducedform
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_dREDUCEDFORM');
                maxrun_dREDUCEDFORM = size(STO_si_dREDUCEDFORM,3);
            else
                maxrun_dREDUCEDFORM = 0;
            end
            if ~options_MC.no_identification_moments
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_dMOMENTS');
                maxrun_dMOMENTS = size(STO_si_dMOMENTS,3);
            else
                maxrun_dMOMENTS = 0;
            end
            if ~options_MC.no_identification_spectrum
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_dSPECTRUM');
                maxrun_dSPECTRUM = size(STO_dSPECTRUM,3);
            else
                maxrun_dSPECTRUM = 0;
            end
            if ~options_MC.no_identification_minimal
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_dMINIMAL');
                maxrun_dMINIMAL = size(STO_dMINIMAL,3);
            else
                maxrun_dMINIMAL = 0;
            end
            for irun=1:max([maxrun_dDYNAMIC, maxrun_dREDUCEDFORM, maxrun_dMOMENTS, maxrun_dSPECTRUM, maxrun_dMINIMAL])
                iter=iter+1;
                % note that this is not the same si_dDYNAMICnorm as computed in identification_analysis
                % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                si_dDYNAMICnorm(iter,:) = vnorm(STO_si_dDYNAMIC(:,:,irun)./repmat(normalize_STO_DYNAMIC,1,totparam_nbr-(stderrparam_nbr+corrparam_nbr))).*normaliz1((stderrparam_nbr+corrparam_nbr)+1:end);
                if ~options_MC.no_identification_reducedform && ~isempty(STO_si_dREDUCEDFORM)
                    % note that this is not the same si_dREDUCEDFORMnorm as computed in identification_analysis
                    % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                    si_dREDUCEDFORMnorm(iter,:) = vnorm(STO_si_dREDUCEDFORM(:,:,irun)./repmat(normalize_STO_REDUCEDFORM,1,totparam_nbr)).*normaliz1;
                end
                if ~options_MC.no_identification_moments && ~isempty(STO_si_dMOMENTS)
                    % note that this is not the same si_dMOMENTSnorm as computed in identification_analysis
                    % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                    si_dMOMENTSnorm(iter,:) = vnorm(STO_si_dMOMENTS(:,:,irun)./repmat(normalize_STO_MOMENTS,1,totparam_nbr)).*normaliz1;
                end
                if ~options_MC.no_identification_spectrum && ~isempty(STO_dSPECTRUM)
                    % note that this is not the same dSPECTRUMnorm as computed in identification_analysis
                    dSPECTRUMnorm(iter,:) = vnorm(STO_dSPECTRUM(:,:,irun)); %not yet used
                end
                if ~options_MC.no_identification_minimal && ~isempty(STO_dMINIMAL)
                    % note that this is not the same dMINIMALnorm as computed in identification_analysis
                    dMINIMALnorm(iter,:) = vnorm(STO_dMINIMAL(:,:,irun)); %not yet used
                end
            end
        end
        IDE_DYNAMIC.si_dDYNAMICnorm = si_dDYNAMICnorm;
        save([IdentifDirectoryName '/' fname '_identif.mat'], 'pdraws', 'IDE_DYNAMIC','STO_DYNAMIC','-append');
        if ~options_MC.no_identification_reducedform
            IDE_REDUCEDFORM.si_dREDUCEDFORMnorm = si_dREDUCEDFORMnorm;
            save([IdentifDirectoryName '/' fname '_identif.mat'], 'IDE_REDUCEDFORM', 'STO_REDUCEDFORM','-append');
        end
        if ~options_MC.no_identification_moments
            IDE_MOMENTS.si_dMOMENTSnorm = si_dMOMENTSnorm;
            save([IdentifDirectoryName '/' fname '_identif.mat'], 'IDE_MOMENTS', 'STO_MOMENTS','-append');
        end

    end

else
    %% load previous analysis
    load([IdentifDirectoryName '/' fname '_identif']);
    parameters                  = store_options_ident.parameter_set;
    options_ident.parameter_set = parameters;
    options_ident.prior_mc      = size(pdraws,1);
    SampleSize                  = options_ident.prior_mc;
    options_.options_ident      = options_ident;
end

%% if dynare_identification is called as it own function (not through identification command) and if we load files
if nargout>3 && iload
    filnam = dir([IdentifDirectoryName '/' fname '_identif_*.mat']);
    STO_si_dDYNAMIC     = [];
    STO_si_dREDUCEDFORM = [];
    STO_si_dMOMENTS     = [];
    STO_dSPECTRUM       = [];
    STO_dMINIMAL        = [];
    for j=1:length(filnam)
        load([IdentifDirectoryName '/' fname '_identif_',int2str(j),'.mat']);
        STO_si_dDYNAMIC = cat(3,STO_si_dDYNAMIC, STO_si_dDYNAMIC(:,abs(iload),:));
        if ~options_ident.no_identification_reducedform
            STO_si_dREDUCEDFORM = cat(3,STO_si_dREDUCEDFORM, STO_si_dREDUCEDFORM(:,abs(iload),:));
        end
        if ~options_ident.no_identification_moments
            STO_si_dMOMENTS = cat(3,STO_si_dMOMENTS, STO_si_dMOMENTS(:,abs(iload),:));
        end
        if ~options_ident.no_identification_spectrum
            STO_dSPECTRUM = cat(3,STO_dSPECTRUM, STO_dSPECTRUM(:,abs(iload),:));
        end
        if ~options_ident.no_identification_minimal
            STO_dMINIMAL = cat(3,STO_dMINIMAL, STO_dMINIMAL(:,abs(iload),:));
        end
    end
end

if iload
    %if previous analysis is loaded
    fprintf(['Testing %s\n',parameters]);
    disp_identification(ide_hess_point.params, ide_reducedform_point, ide_moments_point, ide_spectrum_point, ide_minimal_point, name, options_ident);
    if ~options_.nograph && ~error_indicator_point.identification_strength && ~error_indicator_point.identification_moments
        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
        plot_identification(ide_hess_point.params, ide_moments_point, ide_hess_point, ide_reducedform_point, ide_dynamic_point, options_ident.advanced, parameters, name, IdentifDirectoryName, [], name_tex);
    end
end

%displaying and plotting of results for MC sample
if SampleSize > 1
    fprintf('\nTesting MC sample\n');
    %print results to console but make sure advanced=0
    advanced0 = options_ident.advanced;
    options_ident.advanced = 0;
    disp_identification(pdraws, IDE_REDUCEDFORM, IDE_MOMENTS, IDE_SPECTRUM, IDE_MINIMAL, name, options_ident);
    options_ident.advanced = advanced0; % reset advanced setting
    if ~options_.nograph && isfield(ide_hess_point,'ide_strength_dMOMENTS')
        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
        plot_identification(pdraws, IDE_MOMENTS, ide_hess_point, IDE_REDUCEDFORM, IDE_DYNAMIC, options_ident.advanced, 'MC sample ', name, IdentifDirectoryName, [], name_tex);
    end
    %advanced display and plots for MC Sample, i.e. look at draws with highest/lowest condition number
    if options_ident.advanced
        jcrit = find(IDE_MOMENTS.ino);
        if length(jcrit) < SampleSize
            if isempty(jcrit)
                % Make sure there is no overflow of plots produced (these are saved to the disk anyways)
                store_nodisplay = options_.nodisplay;
                options_.nodisplay = 1;
                % HIGHEST CONDITION NUMBER
                [~, jmax] = max(IDE_MOMENTS.cond);                
                tittxt = 'Draw with HIGHEST condition number';
                fprintf('\nTesting %s.\n',tittxt);
                if ~iload
                    options_ident.tittxt = tittxt; %title text for graphs and figures
                    [ide_moments_max, ide_spectrum_max, ide_minimal_max, ide_hess_max, ide_reducedform_max, ide_dynamic_max, derivatives_info_max, info_max, error_indicator_max] = ...
                        identification_analysis(pdraws(jmax,:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end initializes some persistent variables
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_max', 'ide_moments_max', 'ide_spectrum_max', 'ide_minimal_max','ide_reducedform_max', 'ide_dynamic_max', 'jmax', '-append');
                end
                advanced0 = options_ident.advanced; options_ident.advanced = 1; % make sure advanced setting is on
                disp_identification(pdraws(jmax,:), ide_reducedform_max, ide_moments_max, ide_spectrum_max, ide_minimal_max, name, options_ident);
                options_ident.advanced = advanced0; %reset advanced setting
                if ~options_.nograph && ~error_indicator_max.identification_strength && ~error_indicator_max.identification_moments
                    % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                    plot_identification(pdraws(jmax,:), ide_moments_max, ide_hess_max, ide_reducedform_max, ide_dynamic_max, 1, tittxt, name, IdentifDirectoryName, tittxt, name_tex);
                end

                % SMALLEST condition number
                [~, jmin] = min(IDE_MOMENTS.cond);
                tittxt = 'Draw with SMALLEST condition number';
                fprintf('Testing %s.\n',tittxt);
                if ~iload
                    options_ident.tittxt = tittxt; %title text for graphs and figures
                    [ide_moments_min, ide_spectrum_min, ide_minimal_min, ide_hess_min, ide_reducedform_min, ide_dynamic_min, derivatives_info_min, info_min, error_indicator_min] = ...
                        identification_analysis(pdraws(jmin,:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end initializes persistent variables
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_min', 'ide_moments_min','ide_spectrum_min','ide_minimal_min','ide_reducedform_min', 'ide_dynamic_min', 'jmin', '-append');
                end
                advanced0 = options_ident.advanced; options_ident.advanced = 1; % make sure advanced setting is on
                disp_identification(pdraws(jmin,:), ide_reducedform_min, ide_moments_min, ide_spectrum_min, ide_minimal_min, name, options_ident);
                options_ident.advanced = advanced0; %reset advanced setting
                if ~options_.nograph && ~error_indicator_min.identification_strength && ~error_indicator_min.identification_moments
                    % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                    plot_identification(pdraws(jmin,:),ide_moments_min,ide_hess_min,ide_reducedform_min,ide_dynamic_min,1,tittxt,name,IdentifDirectoryName,tittxt,name_tex);
                end
                % reset nodisplay option
                options_.nodisplay = store_nodisplay;
            else
                % Make sure there is no overflow of plots produced (these are saved to the disk anyways)
                store_nodisplay = options_.nodisplay;
                options_.nodisplay = 1;
                for j=1:length(jcrit)
                    tittxt = ['Rank deficient draw nr. ',int2str(j)];
                    fprintf('\nTesting %s.\n',tittxt);
                    if ~iload
                        options_ident.tittxt = tittxt; %title text for graphs and figures
                        [ide_moments_(j), ide_spectrum_(j), ide_minimal_(j), ide_hess_(j), ide_reducedform_(j), ide_dynamic_(j), derivatives_info_(j), info_resolve, error_indicator_j] = ...
                            identification_analysis(pdraws(jcrit(j),:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1);
                    end
                    advanced0 = options_ident.advanced; options_ident.advanced = 1; %make sure advanced setting is on
                    disp_identification(pdraws(jcrit(j),:), ide_reducedform_(j), ide_moments_(j), ide_spectrum_(j), ide_minimal_(j), name, options_ident);
                    options_ident.advanced = advanced0; % reset advanced
                    if ~options_.nograph && ~error_indicator_j.identification_strength && ~error_indicator_j.identification_moments
                        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                        plot_identification(pdraws(jcrit(j),:), ide_moments_(j), ide_hess_(j), ide_reducedform_(j), ide_dynamic_(j), 1, tittxt, name, IdentifDirectoryName, tittxt, name_tex);
                    end
                end
                if ~iload
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_', 'ide_moments_', 'ide_reducedform_', 'ide_dynamic_', 'ide_spectrum_', 'ide_minimal_', 'jcrit', '-append');
                end
                % reset nodisplay option
                options_.nodisplay = store_nodisplay;
            end
        end
    end
end

%reset warning state
warning_config;

fprintf('\n==== Identification analysis completed ====\n\n')

options_ = store_options_; %restore options set
