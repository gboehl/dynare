function [pdraws, STO_TAU, STO_MOMENTS, STO_LRE, STO_si_dLRE, STO_si_dTAU, STO_si_J, STO_G, STO_D] = dynare_identification(options_ident, pdraws0)
%function [pdraws, STO_TAU, STO_MOMENTS, STO_LRE, STO_dLRE, STO_dTAU, STO_J, STO_G, STO_D] = dynare_identification(options_ident, pdraws0)
% -------------------------------------------------------------------------
% This function is called, when the user specifies identification(...); in
% the mod file. It prepares all identification analysis, i.e.
% (1) sets options, local/persistent/global variables for a new identification 
%     analysis either for a single point or MC Sample and displays and plots the results
% or
% (2) loads, displays and plots a previously saved identification analysis
% =========================================================================
% INPUTS
%    * options_ident    [structure] identification options
%    * pdraws0          [SampleSize by totparam_nbr] optional: matrix of MC sample of model parameters
% -------------------------------------------------------------------------
% OUTPUTS
% Note: This function does not output the arguments to the workspace if only called by
%       "identification" in the mod file, but saves results to the folder identification.
%       One can, however, just use
%       [pdraws, STO_TAU, STO_MOMENTS, STO_LRE, STO_dLRE, STO_dTAU, STO_J, STO_G, STO_D] = dynare_identification(options_ident, pdraws0)
%       in the mod file to get the results directly in the workspace
%    * pdraws           [matrix] MC sample of model params used
%    * STO_TAU,         [matrix] MC sample of entries in the model solution (stacked vertically)
%    * STO_MOMENTS,     [matrix] MC sample of entries in the moments (stacked vertically)
%    * STO_LRE,         [matrix] MC sample of entries in LRE model (stacked vertically)
%    * STO_dLRE,        [matrix] MC sample of derivatives of the Jacobian (dLRE)
%    * STO_dTAU,        [matrix] MC sample of derivatives of the model solution and steady state (dTAU)
%    * STO_J            [matrix] MC sample of Iskrev (2010)'s J matrix
%    * STO_G            [matrix] MC sample of Qu and Tkachenko (2012)'s G matrix
%    * STO_D            [matrix] MC sample of Komunjer and Ng (2011)'s D matrix
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
% Copyright (C) 2010-2019 Dynare Team
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
% =========================================================================

global M_ options_ oo_ bayestopt_ estim_params_

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
end

%% Set all options and create objects
options_ident = set_default_option(options_ident,'gsa_sample_file',0);
    % 0: do not use sample file.
    % 1: triggers gsa prior sample.
    % 2: triggers gsa Monte-Carlo sample (i.e. loads a sample corresponding to pprior=0 and ppost=0 in dynare_sensitivity options).
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
    % 1: use autocorrelations in Iskrev (2010)'s J criteria
    % 0: use autocovariances in Iskrev (2010)'s J criteria
options_ident = set_default_option(options_ident,'ar',1);
    % number of lags to consider for autocovariances/autocorrelations in Iskrev (2010)'s J criteria
options_ident = set_default_option(options_ident,'prior_mc',1);
    % size of Monte-Carlo sample of parameter draws
options_ident = set_default_option(options_ident,'prior_range',0);
    % 1: sample uniformly from prior ranges implied by the prior specifications (overwrites prior shape when prior_mc > 1)
    % 0: sample from specified prior distributions (when prior_mc > 1)
options_ident = set_default_option(options_ident,'periods',300);
    % length of stochastic simulation to compute simulated moment uncertainty, when analytic Hessian is not available
options_ident = set_default_option(options_ident,'replic',100);
    % number of replicas to compute simulated moment uncertainty, when analytic Hessian is not available
options_ident = set_default_option(options_ident,'advanced',0);
    % 1: show a more detailed analysis based on reduced-form solution and Jacobian of dynamic model (LRE). Further, performs a brute force
    %    search of the groups of parameters best reproducing the behavior of each single parameter of Iskrev (2010)'s J.
options_ident = set_default_option(options_ident,'normalize_jacobians',1);
    % 1: normalize Jacobians by rescaling each row by its largest element in absolute value
options_ident = set_default_option(options_ident,'grid_nbr',5000);
    % number of grid points in [-pi;pi] to approximate the integral to compute Qu and Tkachenko (2012)'s G criteria
    % note that grid_nbr needs to be even and actually we use (grid_nbr+1) points, as we add the 0 frequency and use symmetry, i.e. grid_nbr/2
    % negative as well as grid_nbr/2 positive values to speed up the compuations
    if mod(options_ident.grid_nbr,2) ~= 0
        options_ident.grid_nbr = options_ident.grid_nbr+1;
        if mod(options_ident.grid_nbr,2) ~= 0
            error('IDENTIFICATION: You need to set an even value for ''grid_nbr''');
        end
    end
options_ident = set_default_option(options_ident,'tol_rank',1.e-10);
    % tolerance level used for rank computations
options_ident = set_default_option(options_ident,'tol_deriv',1.e-8);
    % tolerance level for selecting columns of non-zero derivatives
options_ident = set_default_option(options_ident,'tol_sv',1.e-3);
    % tolerance level for selecting non-zero singular values in identification_checks.m
    
%check whether to compute identification strength based on information matrix
if ~isfield(options_ident,'no_identification_strength')
    options_ident.no_identification_strength = 0;
else
    options_ident.no_identification_strength = 1;
end
%check whether to compute and display identification criteria based on steady state and reduced-form solution
if ~isfield(options_ident,'no_identification_reducedform')
    options_ident.no_identification_reducedform = 0;
else
    options_ident.no_identification_reducedform = 1;
end
%check whether to compute and display identification criteria based on Iskrev (2010)'s J, i.e. derivative of first two moments
if ~isfield(options_ident,'no_identification_moments')
    options_ident.no_identification_moments = 0;
else
    options_ident.no_identification_moments = 1;
end
%check whether to compute and display identification criteria based on Komunjer and Ng (2011)'s D, i.e. derivative of first moment, minimal state space system and observational equivalent spectral density transformation
if ~isfield(options_ident,'no_identification_minimal')
    options_ident.no_identification_minimal = 0;
else
     options_ident.no_identification_minimal = 1;
end
%Check whether to compute and display identification criteria based on Qu and Tkachenko (2012)'s G, i.e. Gram matrix of derivatives of first moment plus outer product of derivatives of spectral density
if ~isfield(options_ident,'no_identification_spectrum')
    options_ident.no_identification_spectrum = 0;
else
    options_ident.no_identification_spectrum = 1;
end

%overwrite setting, as identification strength and advanced need criteria based on both reducedform and moments
if (isfield(options_ident,'no_identification_strength') &&  options_ident.no_identification_strength == 0) || (options_ident.advanced == 1)
    options_ident.no_identification_reducedform = 0;
    options_ident.no_identification_moments = 0;
end

%overwrite setting, as dynare_sensitivity does not make use of spectrum and minimal system (yet)
if isfield(options_,'opt_gsa') && isfield(options_.opt_gsa,'identification') && options_.opt_gsa.identification == 1
    options_ident.no_identification_minimal = 1;
    options_ident.no_identification_spectrum = 1;
end

%Deal with non-stationary cases
if isfield(options_ident,'diffuse_filter') %set lik_init and options_
    options_ident.lik_init=3; %overwrites any other lik_init set
    options_.diffuse_filter=options_ident.diffuse_filter; %make options_ inherit diffuse filter; will overwrite any conflicting lik_init in dynare_estimation_init
else
    if options_.diffuse_filter==1 %warning if estimation with diffuse filter was done, but not passed
        fprintf('WARNING IDENTIFICATION: Previously the diffuse_filter option was used, but it was not passed to the identification command. This may result in problems if your model contains unit roots.\n');
    end
    if isfield(options_ident,'lik_init')
        options_.lik_init=options_ident.lik_init; %make options_ inherit lik_init
        if options_ident.lik_init==3 %user specified diffuse filter using the lik_init option
            options_ident.analytic_derivation=0; %diffuse filter not compatible with analytic derivation
            options_.analytic_derivation=0; %diffuse filter not compatible with analytic derivation
        end
    end
end
options_ident = set_default_option(options_ident,'lik_init',1);
    % Type of initialization of Kalman filter:
    % 1: stationary models: initial matrix of variance of error of forecast is set equal to the unconditional variance of the state variables
    % 2: nonstationary models: wide prior is used with an initial matrix of variance of the error of forecast diagonal with 10 on the diagonal (follows the suggestion of Harvey and Phillips(1979))
    % 3: nonstationary models: use a diffuse filter (use rather the diffuse_filter option)
    % 4: filter is initialized with the fixed point of the Riccati equation
    % 5:  i) option 2 for non-stationary elements by setting their initial variance in the forecast error matrix to 10 on the diagonal and all co-variances to 0 and
    %    ii) option 1 for the stationary elements
options_ident = set_default_option(options_ident,'analytic_derivation',1);
    % 1: analytic derivation of gradient and hessian of likelihood in dsge_likelihood.m, only works for stationary models, i.e. kalman_algo<3

%overwrite values in options_, note that options_ is restored at the end of the function
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
else
    prior_exist = 1;
    parameters = options_ident.parameter_set;
end

% overwrite settings in options_ and prepare to call dynare_estimation_init
options_.order = 1; % Identification does not support order>1 (yet)
options_.ar = options_ident.ar;
options_.prior_mc = options_ident.prior_mc;
options_.Schur_vec_tol = 1.e-8;
options_.nomoments = 0;
options_ = set_default_option(options_,'analytic_derivation',1); %if option was not already set
    % 1: analytic derivation of gradient and hessian of likelihood in dsge_likelihood.m, only works for stationary models, i.e. kalman_algo<3
options_ = set_default_option(options_,'datafile','');
options_.mode_compute = 0;
options_.plot_priors = 0;
options_.smoother = 1;
[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = dynare_estimation_init(M_.endo_names, fname, 1, M_, options_, oo_, estim_params_, bayestopt_);
%outputting dataset_ is needed for Octave

% set method to compute identification Jacobians (kronflag). Default:0
options_ident = set_default_option(options_ident,'analytic_derivation_mode', options_.analytic_derivation_mode); % if not set by user, inherit default global one
    %  0: efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2011)
    %  1: kronecker products method to compute analytical derivatives as in Iskrev (2010)
    % -1: numerical two-sided finite difference method to compute numerical derivatives of all Jacobians using function identification_numerical_objective.m (previously thet2tau.m)
    % -2: numerical two-sided finite difference method to compute numerically dYss, dg1, d2Yss and d2g1, the Jacobians are then computed analytically as in options 0
options_.analytic_derivation_mode = options_ident.analytic_derivation_mode; %overwrite setting

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
    modparam_nbr = estim_params_.np; %number of model parameters as declared in estimated_params
    stderrparam_nbr = estim_params_.nvx; % nvx is number of stderr parameters
    corrparam_nbr = estim_params_.ncx; % ncx is number of corr parameters
    if estim_params_.nvn || estim_params_.ncn %nvn is number of stderr parameters and ncn is number of corr parameters of measurement innovations as declared in estimated_params
        error('Identification does not (yet) support measurement errors. Instead, define them explicitly in measurement equations in model definition.')
    end
    name = cell(totparam_nbr,1); %initialize cell for parameter names
    name_tex = cell(totparam_nbr,1);
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
    indpmodel = 1:M_.param_nbr; %all model  parameters
    indpstderr = 1:M_.exo_nbr;  %all stderr parameters
    indpcorr = [];              %no  corr   parameters
    stderrparam_nbr = M_.exo_nbr;
    corrparam_nbr = 0;
    modparam_nbr = M_.param_nbr;
    totparam_nbr = modparam_nbr+stderrparam_nbr;
    name = cellfun(@(x) horzcat('SE_', x), M_.exo_names, 'UniformOutput', false); %names for stderr parameters
    name = vertcat(name, M_.param_names);
    name_tex = cellfun(@(x) horzcat('$ SE_{', x, '} $'), M_.exo_names, 'UniformOutput', false);
    name_tex = vertcat(name_tex, M_.param_names_tex);
    if ~isequal(M_.H,0)
        fprintf('\ndynare_identification:: Identification does not support measurement errors (yet) and will ignore them in the following. To test their identifiability, instead define them explicitly in measurement equations in the model definition.\n')
    end
end
options_ident.name_tex = name_tex;

skipline()
disp('==== Identification analysis ====')
skipline()
if totparam_nbr < 2
    options_ident.advanced = 0;
    disp('There is only one parameter to study for identitification. The advanced option is re-set to 0.')
    skipline()
end

% set options_ident dependent ot totparam_nbr
options_ident = set_default_option(options_ident,'max_dim_cova_group',min([2,totparam_nbr-1]));
options_ident.max_dim_cova_group = min([options_ident.max_dim_cova_group,totparam_nbr-1]);
    % In brute force search (ident_bruteforce.m) when advanced=1 this option sets the maximum dimension of groups of parameters that best reproduce the behavior of each single model parameter

options_ident = set_default_option(options_ident,'checks_via_subsets',0); %[ONLY FOR DEBUGGING]
    % 1: uses identification_checks_via_subsets.m to compute problematic parameter combinations 
    % 0: uses identification_checks.m to compute problematic parameter combinations [default]
options_ident = set_default_option(options_ident,'max_dim_subsets_groups',min([4,totparam_nbr-1])); %[ONLY FOR DEBUGGING]    
    % In identification_checks_via_subsets.m, when checks_via_subsets=1, this
    % option sets the maximum dimension of groups of parameters for which
    % the corresponding rank criteria is checked 

options_.options_ident = options_ident; %store identification options into global options
MaxNumberOfBytes = options_.MaxNumberOfBytes; %threshold when to save from memory to files
store_options_ident = options_ident;

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
            disp('Testing ML Starting value')
        else
            % use user-defined option
            switch parameters
              case 'calibration'
                parameters_TeX = 'Calibration';
                disp('Testing calibration')
                params(1,:) = get_all_parameters(estim_params_,M_);
              case 'posterior_mode'
                parameters_TeX = 'Posterior mode';
                disp('Testing posterior mode')
                params(1,:) = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
              case 'posterior_mean'
                parameters_TeX = 'Posterior mean';
                disp('Testing posterior mean')
                params(1,:) = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
              case 'posterior_median'
                parameters_TeX = 'Posterior median';
                disp('Testing posterior median')
                params(1,:) = get_posterior_parameters('median',M_,estim_params_,oo_,options_);
              case 'prior_mode'
                parameters_TeX = 'Prior mode';
                disp('Testing prior mode')
                params(1,:) = bayestopt_.p5(:);
              case 'prior_mean'
                parameters_TeX = 'Prior mean';
                disp('Testing prior mean')
                params(1,:) = bayestopt_.p1;
              otherwise
                disp('The option parameter_set has to be equal to:')
                disp('                   ''calibration'', ')
                disp('                   ''posterior_mode'', ')
                disp('                   ''posterior_mean'', ')
                disp('                   ''posterior_median'', ')
                disp('                   ''prior_mode'' or')
                disp('                   ''prior_mean''.')
                error('IDENTIFICATION: The option ''parameter_set'' has and invalid value');
            end
        end
    else
        % no estimated_params block is available, all stderr and model parameters, but no corr parameters are chosen
        params = [sqrt(diag(M_.Sigma_e))', M_.params']; % use current values
        parameters = 'Current_params';
        parameters_TeX = 'Current parameter values';
        disp('Testing all current stderr and model parameter values')
    end
    options_ident.tittxt = parameters; %title text for graphs and figures
    % perform identification analysis for single point
    [ide_moments_point, ide_spectrum_point, ide_minimal_point, ide_hess_point, ide_reducedform_point, ide_lre_point, derivatives_info_point, info, options_ident] = ...
        identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end implies initialization of persistent variables
    if info(1)~=0
        % there are errors in the solution algorithm
        skipline()
        disp('----------- ')
        disp('Parameter error:')
        disp(['The model does not solve for ', parameters, ' with error code info = ', int2str(info(1))]),
        skipline()
        if info(1)==1
            disp('info==1 %! The model doesn''t determine the current variables uniquely.')
        elseif info(1)==2
            disp('info==2 %! MJDGGES returned an error code.')
        elseif info(1)==3
            disp('info==3 %! Blanchard & Kahn conditions are not satisfied: no stable equilibrium. ')
        elseif info(1)==4
            disp('info==4 %! Blanchard & Kahn conditions are not satisfied: indeterminacy. ')
        elseif info(1)==5
            disp('info==5 %! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure. ')
        elseif info(1)==6
            disp('info==6 %! The jacobian evaluated at the deterministic steady state is complex.')
        elseif info(1)==19
            disp('info==19 %! The steadystate routine has thrown an exception (inconsistent deep parameters). ')
        elseif info(1)==20
            disp('info==20 %! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations). ')
        elseif info(1)==21
            disp('info==21 %! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.')
        elseif info(1)==22
            disp('info==22 %! The steady has NaNs. ')
        elseif info(1)==23
            disp('info==23 %! M_.params has been updated in the steadystate routine and has complex valued scalars. ')
        elseif info(1)==24
            disp('info==24 %! M_.params has been updated in the steadystate routine and has some NaNs. ')
        elseif info(1)==30
            disp('info==30 %! Ergodic variance can''t be computed. ')
        end
        disp('----------- ')
        skipline()
        if any(bayestopt_.pshape)
            % if there are errors in the solution algorithm, try to sample a different point from the prior
            disp('Try sampling up to 50 parameter sets from the prior.')
            kk=0;
            while kk<50 && info(1)
                kk=kk+1;
                params = prior_draw();
                options_ident.tittxt = 'Random_prior_params'; %title text for graphs and figures
                % perform identification analysis
                [ide_moments_point, ide_spectrum_point, ide_minimal_point, ide_hess_point, ide_reducedform_point, ide_lre_point, derivatives_info, info, options_ident] = ...
                    identification_analysis(params,indpmodel,indpstderr,indpcorr,options_ident,dataset_info, prior_exist, 1);
            end
        end
        if info(1)
            skipline()
            disp('----------- ')
            disp('Identification stopped:')
            if any(bayestopt_.pshape)
                disp('The model did not solve for any of 50 attempts of random samples from the prior')
            end
            disp('----------- ')
            skipline()
            return
        else
            % found a (random) point that solves the model
            disp('Found a random draw from the priors that solves the model.')
            disp(params)
            disp('Identification now continues for this draw.');
            parameters = 'Random_prior_params';
            parameters_TeX = 'Random prior parameter draw';
        end
    end
    ide_hess_point.params = params;
    % save all output into identification folder
    save([IdentifDirectoryName '/' fname                '_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point', 'ide_hess_point', 'ide_reducedform_point', 'ide_lre_point','store_options_ident');
    save([IdentifDirectoryName '/' fname '_' parameters '_identif.mat'], 'ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point', 'ide_hess_point', 'ide_reducedform_point', 'ide_lre_point','store_options_ident');
    % display results of identification analysis
    disp_identification(params, ide_reducedform_point, ide_moments_point, ide_spectrum_point, ide_minimal_point, name, options_ident);
    if ~options_ident.no_identification_strength && ~options_.nograph
        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
        plot_identification(params, ide_moments_point, ide_hess_point, ide_reducedform_point, ide_lre_point, options_ident.advanced, parameters, name, IdentifDirectoryName, parameters_TeX, name_tex);
    end
    
    if SampleSize > 1
        % initializations for Monte Carlo Analysis
        skipline()
        disp('Monte Carlo Testing')
        h = dyn_waitbar(0,'Monte Carlo identification checks ...');
        iteration  = 0; % initialize counter for admissable draws
        run_index  = 0; % initialize counter for admissable draws after saving previous draws to file(s)
        file_index = 0; % initialize counter for files (if MaxNumberOfBytes is reached, we store results in files)
        options_MC = options_ident; %store options structure for Monte Carlo analysis
        options_MC.advanced = 0;    %do not run advanced checking in a Monte Carlo analysis        
        options_ident.checks_via_subsets = 0; % for Monte Carlo analysis currently only identification_checks and not identification_checks_via_subsets is supported
    else
        iteration = 1; % iteration equals SampleSize and we are finished
        pdraws = [];   % to have output object otherwise map_ident may crash
    end
    while iteration < SampleSize
        if external_sample
            params = pdraws0(iteration+1,:); % loaded draws
        else
            params = prior_draw(); % new random draw from prior
        end
        options_ident.tittxt = []; % clear title text for graphs and figures
        % run identification analysis
        [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_lre, ide_derivatives_info, info, options_MC] = ...
            identification_analysis(params, indpmodel, indpstderr, indpcorr, options_MC, dataset_info, prior_exist, 0); % the 0 implies that we do not initialize persistent variables anymore
        
        if iteration==0 && info(1)==0 % preallocate storage in the first admissable run
            delete([IdentifDirectoryName '/' fname '_identif_*.mat']) % delete previously saved results
            MAX_RUNS_BEFORE_SAVE_TO_FILE = min(SampleSize,ceil(MaxNumberOfBytes/(size(ide_reducedform.si_dTAU,1)*totparam_nbr)/8)); % set how many runs can be stored before we save to files
            pdraws = zeros(SampleSize,totparam_nbr); % preallocate storage for draws in each row
            
            % preallocate storage for linear rational expectations model
            STO_si_dLRE        = zeros([size(ide_lre.si_dLRE,1),modparam_nbr,MAX_RUNS_BEFORE_SAVE_TO_FILE]);
            STO_LRE            = zeros(size(ide_lre.LRE,1),SampleSize);
            IDE_LRE.ind_dLRE   = ide_lre.ind_dLRE;
            IDE_LRE.in0        = zeros(SampleSize,modparam_nbr);
            IDE_LRE.ind0       = zeros(SampleSize,modparam_nbr);
            IDE_LRE.jweak      = zeros(SampleSize,modparam_nbr);
            IDE_LRE.jweak_pair = zeros(SampleSize,modparam_nbr*(modparam_nbr+1)/2);
            IDE_LRE.cond       = zeros(SampleSize,1);
            IDE_LRE.Mco        = zeros(SampleSize,modparam_nbr);
            
            % preallocate storage for reduced form
            if ~options_MC.no_identification_reducedform
                STO_si_dTAU                = zeros([size(ide_reducedform.si_dTAU,1),totparam_nbr,MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                STO_TAU                    = zeros(size(ide_reducedform.TAU,1),SampleSize);
                IDE_REDUCEDFORM.ind_dTAU   = ide_reducedform.ind_dTAU;
                IDE_REDUCEDFORM.in0        = zeros(SampleSize,1);
                IDE_REDUCEDFORM.ind0       = zeros(SampleSize,totparam_nbr);
                IDE_REDUCEDFORM.jweak      = zeros(SampleSize,totparam_nbr);
                IDE_REDUCEDFORM.jweak_pair = zeros(SampleSize,totparam_nbr*(totparam_nbr+1)/2);
                IDE_REDUCEDFORM.cond       = zeros(SampleSize,1);
                IDE_REDUCEDFORM.Mco        = zeros(SampleSize,totparam_nbr);
            else
                IDE_REDUCEDFORM = {};
            end
            
            % preallocate storage for moments
            if ~options_MC.no_identification_moments
                STO_si_J               = zeros([size(ide_moments.si_J,1),totparam_nbr,MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                STO_MOMENTS            = zeros(size(ide_moments.MOMENTS,1),SampleSize);
                IDE_MOMENTS.ind_J      = ide_moments.ind_J;
                IDE_MOMENTS.in0        = zeros(SampleSize,1);
                IDE_MOMENTS.ind0       = zeros(SampleSize,totparam_nbr);
                IDE_MOMENTS.jweak      = zeros(SampleSize,totparam_nbr);
                IDE_MOMENTS.jweak_pair = zeros(SampleSize,totparam_nbr*(totparam_nbr+1)/2);
                IDE_MOMENTS.cond       = zeros(SampleSize,1);
                IDE_MOMENTS.Mco        = zeros(SampleSize,totparam_nbr);
                IDE_MOMENTS.S          = zeros(SampleSize,min(8,totparam_nbr));
                IDE_MOMENTS.V          = zeros(SampleSize,totparam_nbr,min(8,totparam_nbr));
            else
                IDE_MOMENTS = {};
            end
            
            % preallocate storage for spectrum
            if ~options_MC.no_identification_spectrum
                STO_G                   = zeros([size(ide_spectrum.G,1),size(ide_spectrum.G,2), MAX_RUNS_BEFORE_SAVE_TO_FILE]);
                IDE_SPECTRUM.ind_G      = ide_spectrum.ind_G;
                IDE_SPECTRUM.in0        = zeros(SampleSize,1);
                IDE_SPECTRUM.ind0       = zeros(SampleSize,totparam_nbr);
                IDE_SPECTRUM.jweak      = zeros(SampleSize,totparam_nbr);
                IDE_SPECTRUM.jweak_pair = zeros(SampleSize,totparam_nbr*(totparam_nbr+1)/2);
                IDE_SPECTRUM.cond       = zeros(SampleSize,1);
                IDE_SPECTRUM.Mco        = zeros(SampleSize,totparam_nbr);
            else
                IDE_SPECTRUM = {};
            end
            
            % preallocate storage for minimal system
            if ~options_MC.no_identification_minimal
                STO_D                  = zeros([size(ide_minimal.D,1),size(ide_minimal.D,2), MAX_RUNS_BEFORE_SAVE_TO_FILE]);                
                IDE_MINIMAL.ind_D      = ide_minimal.ind_D;
                IDE_MINIMAL.in0        = zeros(SampleSize,1);
                IDE_MINIMAL.ind0       = zeros(SampleSize,totparam_nbr);
                IDE_MINIMAL.jweak      = zeros(SampleSize,totparam_nbr);
                IDE_MINIMAL.jweak_pair = zeros(SampleSize,totparam_nbr*(totparam_nbr+1)/2);
                IDE_MINIMAL.cond       = zeros(SampleSize,1);
                IDE_MINIMAL.Mco        = zeros(SampleSize,totparam_nbr);
            else
                IDE_MINIMAL = {};
            end
        end
        
        if info(1)==0 % if admissable draw
            iteration = iteration + 1; %increase total index of admissable draws
            run_index = run_index + 1; %increase index of admissable draws after saving to files
            pdraws(iteration,:) = params; % store draw
            
            % store results for linear rational expectations model
            STO_LRE(:,iteration)            = ide_lre.LRE;
            STO_si_dLRE(:,:,run_index)      = ide_lre.si_dLRE;
            IDE_LRE.cond(iteration,1)       = ide_lre.cond;
            IDE_LRE.ino(iteration,1)        = ide_lre.ino;
            IDE_LRE.ind0(iteration,:)       = ide_lre.ind0;
            IDE_LRE.jweak(iteration,:)      = ide_lre.jweak;
            IDE_LRE.jweak_pair(iteration,:) = ide_lre.jweak_pair;
            IDE_LRE.Mco(iteration,:)        = ide_lre.Mco;            
            
            % store results for reduced form solution
            if ~options_MC.no_identification_reducedform
                STO_TAU(:,iteration)                    = ide_reducedform.TAU;
                STO_si_dTAU(:,:,run_index)              = ide_reducedform.si_dTAU;
                IDE_REDUCEDFORM.cond(iteration,1)       = ide_reducedform.cond;
                IDE_REDUCEDFORM.ino(iteration,1)        = ide_reducedform.ino;
                IDE_REDUCEDFORM.ind0(iteration,:)       = ide_reducedform.ind0;
                IDE_REDUCEDFORM.jweak(iteration,:)      = ide_reducedform.jweak;
                IDE_REDUCEDFORM.jweak_pair(iteration,:) = ide_reducedform.jweak_pair;
                IDE_REDUCEDFORM.Mco(iteration,:)        = ide_reducedform.Mco;
            end
            
            % store results for moments
            if ~options_MC.no_identification_moments
                STO_MOMENTS(:,iteration)            = ide_moments.MOMENTS;
                STO_si_J(:,:,run_index)             = ide_moments.si_J;
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
            if ~options_MC.no_identification_spectrum
                STO_G(:,:,run_index)                 = ide_spectrum.G;
                IDE_SPECTRUM.cond(iteration,1)       = ide_spectrum.cond;
                IDE_SPECTRUM.ino(iteration,1)        = ide_spectrum.ino;
                IDE_SPECTRUM.ind0(iteration,:)       = ide_spectrum.ind0;
                IDE_SPECTRUM.jweak(iteration,:)      = ide_spectrum.jweak;
                IDE_SPECTRUM.jweak_pair(iteration,:) = ide_spectrum.jweak_pair;
                IDE_SPECTRUM.Mco(iteration,:)        = ide_spectrum.Mco;                
            end
            
            % store results for minimal system
            if ~options_MC.no_identification_minimal
                STO_D(:,:,run_index)                = ide_minimal.D;
                IDE_MINIMAL.cond(iteration,1)       = ide_minimal.cond;
                IDE_MINIMAL.ino(iteration,1)        = ide_minimal.ino;
                IDE_MINIMAL.ind0(iteration,:)       = ide_minimal.ind0;
                IDE_MINIMAL.jweak(iteration,:)      = ide_minimal.jweak;
                IDE_MINIMAL.jweak_pair(iteration,:) = ide_minimal.jweak_pair;
                IDE_MINIMAL.Mco(iteration,:)        = ide_minimal.Mco;                
            end
            
            % save results to file: either to avoid running into memory issues, i.e. (run_index==MAX_RUNS_BEFORE_SAVE_TO_FILE) or if finished (iteration==SampleSize)
            if run_index==MAX_RUNS_BEFORE_SAVE_TO_FILE || iteration==SampleSize
                file_index = file_index + 1;
                if run_index<MAX_RUNS_BEFORE_SAVE_TO_FILE 
                    %we are finished (iteration == SampleSize), so get rid of additional storage
                    STO_si_dLRE = STO_si_dLRE(:,:,1:run_index);
                    if ~options_MC.no_identification_reducedform
                        STO_si_dTAU = STO_si_dTAU(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_moments
                        STO_si_J = STO_si_J(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_spectrum
                        STO_G = STO_G(:,:,1:run_index);
                    end
                    if ~options_MC.no_identification_minimal
                        STO_D = STO_D(:,:,1:run_index);
                    end
                end                
                save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_dLRE');
                STO_si_dLRE = zeros(size(STO_si_dLRE)); % reset storage
                if ~options_MC.no_identification_reducedform
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_dTAU', '-append');
                    STO_si_dTAU = zeros(size(STO_si_dTAU)); % reset storage
                end
                if ~options_MC.no_identification_moments
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_si_J','-append');
                    STO_si_J = zeros(size(STO_si_J)); % reset storage
                end
                if ~options_MC.no_identification_spectrum
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_G','-append');
                    STO_G = zeros(size(STO_G)); % reset storage
                end
                if ~options_MC.no_identification_minimal
                    save([IdentifDirectoryName '/' fname '_identif_' int2str(file_index) '.mat'], 'STO_D','-append');
                    STO_D = zeros(size(STO_D)); % reset storage
                end
                run_index = 0; % reset index
            end
            if SampleSize > 1
                dyn_waitbar(iteration/SampleSize,h,['MC identification checks ',int2str(iteration),'/',int2str(SampleSize)])
            end
        end
    end

    if SampleSize > 1
        dyn_waitbar_close(h);
        normalize_STO_LRE = std(STO_LRE,0,2);
        if ~options_MC.no_identification_reducedform
            normalize_STO_TAU = std(STO_TAU,0,2);
        end
        if ~options_MC.no_identification_moments
            normalize_STO_MOMENTS = std(STO_MOMENTS,0,2);
        end
        if ~options_MC.no_identification_minimal
            normalize_STO_MINIMAL = 1; %not yet used
        end
        if ~options_MC.no_identification_spectrum
            normalize_STO_SPECTRUM = 1; %not yet used
        end
        normaliz1 = std(pdraws);
        iter = 0;
        for ifile_index = 1:file_index
            load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_dLRE');
            maxrun_dLRE = size(STO_si_dLRE,3);
            if ~options_MC.no_identification_reducedform
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_dTAU');
                maxrun_dTAU = size(STO_si_dTAU,3);
            else
                maxrun_dTAU = 0;
            end
            if ~options_MC.no_identification_moments
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_si_J');
                maxrun_J = size(STO_si_J,3);
            else
                maxrun_J = 0;
            end
            if ~options_MC.no_identification_spectrum
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_G');
                maxrun_G = size(STO_G,3);
            else
                maxrun_G = 0;
            end
            if ~options_MC.no_identification_minimal
                load([IdentifDirectoryName '/' fname '_identif_' int2str(ifile_index) '.mat'], 'STO_D');
                maxrun_D = size(STO_D,3);
            else
                maxrun_D = 0;
            end
            for irun=1:max([maxrun_dLRE, maxrun_dTAU, maxrun_J, maxrun_G, maxrun_D])
                iter=iter+1;
                % note that this is not the same si_dLREnorm as computed in identification_analysis
                % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                si_dLREnorm(iter,:) = vnorm(STO_si_dLRE(:,:,irun)./repmat(normalize_STO_LRE,1,totparam_nbr-(stderrparam_nbr+corrparam_nbr))).*normaliz1((stderrparam_nbr+corrparam_nbr)+1:end);
                if ~options_MC.no_identification_reducedform
                    % note that this is not the same si_dTAUnorm as computed in identification_analysis
                    % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                    si_dTAUnorm(iter,:) = vnorm(STO_si_dTAU(:,:,irun)./repmat(normalize_STO_TAU,1,totparam_nbr)).*normaliz1;
                end
                if ~options_MC.no_identification_moments
                    % note that this is not the same si_Jnorm as computed in identification_analysis
                    % given that we have the MC sample of the Jacobians, we also normalize by the std of the sample of Jacobian entries, to get a fully standardized sensitivity measure
                    si_Jnorm(iter,:) = vnorm(STO_si_J(:,:,irun)./repmat(normalize_STO_MOMENTS,1,totparam_nbr)).*normaliz1;
                end
                if ~options_MC.no_identification_spectrum
                    % note that this is not the same Gnorm as computed in identification_analysis                    
                    Gnorm(iter,:) = vnorm(STO_G(:,:,irun)); %not yet used
                end
                if ~options_MC.no_identification_minimal
                    % note that this is not the same Dnorm as computed in identification_analysis
                    Dnorm(iter,:) = vnorm(STO_D(:,:,irun)); %not yet used
                end
            end
        end
        IDE_LRE.si_dLREnorm = si_dLREnorm;
        save([IdentifDirectoryName '/' fname '_identif.mat'], 'pdraws', 'IDE_LRE','STO_LRE','-append');
        if ~options_MC.no_identification_reducedform
            IDE_REDUCEDFORM.si_dTAUnorm = si_dTAUnorm;
            save([IdentifDirectoryName '/' fname '_identif.mat'], 'IDE_REDUCEDFORM', 'STO_TAU','-append');
        end
        if ~options_MC.no_identification_moments
            IDE_MOMENTS.si_Jnorm = si_Jnorm;
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
    STO_si_dLRE = [];
    STO_si_dTAU=[];    
    STO_si_J = [];
    STO_G = [];
    STO_D = [];    
    for j=1:length(filnam)
        load([IdentifDirectoryName '/' fname '_identif_',int2str(j),'.mat']);
        STO_si_dLRE = cat(3,STO_si_dLRE, STO_si_dLRE(:,abs(iload),:));
        if ~options_ident.no_identification_reducedform
            STO_si_dTAU = cat(3,STO_si_dTAU, STO_si_dTAU(:,abs(iload),:));
        end
        if ~options_ident.no_identification_moments
            STO_si_J = cat(3,STO_si_J, STO_si_J(:,abs(iload),:));
        end
        if ~options_ident.no_identification_spectrum
            STO_G = cat(3,STO_G, STO_G(:,abs(iload),:));
        end
        if ~options_ident.no_identification_minimal
            STO_D = cat(3,STO_D, STO_D(:,abs(iload),:));
        end
    end
end

if iload
    %if previous analysis is loaded
    disp(['Testing ',parameters])
    disp_identification(ide_hess_point.params, ide_reducedform_point, ide_moments_point, ide_spectrum_point, ide_minimal_point, name, options_ident);
    if ~options_.nograph
        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
        plot_identification(ide_hess_point.params, ide_moments_point, ide_hess_point, ide_reducedform_point, ide_lre_point, options_ident.advanced, parameters, name, IdentifDirectoryName, [],name_tex);
    end
end

%displaying and plotting of results for MC sample
if SampleSize > 1    
    fprintf('\n')
    disp('Testing MC sample')
    %print results to console but make sure advanced=0
    advanced0 = options_ident.advanced;
    options_ident.advanced = 0;
    disp_identification(pdraws, IDE_REDUCEDFORM, IDE_MOMENTS, IDE_SPECTRUM, IDE_MINIMAL, name, options_ident);
    options_ident.advanced = advanced0; % reset advanced setting
    if ~options_.nograph
        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
        plot_identification(pdraws, IDE_MOMENTS, ide_hess_point, IDE_REDUCEDFORM, IDE_LRE, options_ident.advanced, 'MC sample ', name, IdentifDirectoryName, [], name_tex);
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
                fprintf('\n')
                tittxt = 'Draw with HIGHEST condition number';
                fprintf('\n')
                disp(['Testing ',tittxt, '.']),
                if ~iload
                    options_ident.tittxt = tittxt; %title text for graphs and figures
                    [ide_moments_max, ide_spectrum_max, ide_minimal_max, ide_hess_max, ide_reducedform_max, ide_lre_max, derivatives_info_max, info_max, options_ident] = ...
                        identification_analysis(pdraws(jmax,:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end initializes some persistent variables
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_max', 'ide_moments_max', 'ide_spectrum_max', 'ide_minimal_max','ide_reducedform_max', 'ide_lre_max', 'jmax', '-append');
                end                
                advanced0 = options_ident.advanced; options_ident.advanced = 1; % make sure advanced setting is on
                disp_identification(pdraws(jmax,:), ide_reducedform_max, ide_moments_max, ide_spectrum_max, ide_minimal_max, name, options_ident);
                options_ident.advanced = advanced0; %reset advanced setting                
                if ~options_.nograph
                    % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                    plot_identification(pdraws(jmax,:), ide_moments_max, ide_hess_max, ide_reducedform_max, ide_lre_max, 1, tittxt, name, IdentifDirectoryName, tittxt, name_tex);
                end
                
                % SMALLEST condition number
                [~, jmin] = min(IDE_MOMENTS.cond);
                fprintf('\n')
                tittxt = 'Draw with SMALLEST condition number';
                fprintf('\n')
                disp(['Testing ',tittxt, '.']),
                if ~iload
                    options_ident.tittxt = tittxt; %title text for graphs and figures
                    [ide_moments_min, ide_spectrum_min, ide_minimal_min, ide_hess_min, ide_reducedform_min, ide_lre_min, derivatives_info_min, info_min, options_ident] = ...
                        identification_analysis(pdraws(jmin,:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1); %the 1 at the end initializes persistent variables
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_min', 'ide_moments_min','ide_spectrum_min','ide_minimal_min','ide_reducedform_min', 'ide_lre_min', 'jmin', '-append');
                end
                advanced0 = options_ident.advanced; options_ident.advanced = 1; % make sure advanced setting is on
                disp_identification(pdraws(jmin,:), ide_reducedform_min, ide_moments_min, ide_spectrum_min, ide_minimal_min, name, options_ident);
                options_ident.advanced = advanced0; %reset advanced setting
                if ~options_.nograph
                    % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                    plot_identification(pdraws(jmin,:),ide_moments_min,ide_hess_min,ide_reducedform_min,ide_lre_min,1,tittxt,name,IdentifDirectoryName,tittxt,name_tex);
                end
                % reset nodisplay option
                options_.nodisplay = store_nodisplay;
            else
                % Make sure there is no overflow of plots produced (these are saved to the disk anyways)
                store_nodisplay = options_.nodisplay;
                options_.nodisplay = 1;
                for j=1:length(jcrit)
                    tittxt = ['Rank deficient draw n. ',int2str(j)];
                    fprintf('\n')
                    disp(['Testing ',tittxt, '.']),
                    if ~iload
                        options_ident.tittxt = tittxt; %title text for graphs and figures
                        [ide_moments_(j), ide_spectrum_(j), ide_minimal_(j), ide_hess_(j), ide_reducedform_(j), ide_lre_(j), derivatives_info_(j), info_resolve, options_ident] = ...
                            identification_analysis(pdraws(jcrit(j),:), indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, 1);
                    end
                    advanced0 = options_ident.advanced; options_ident.advanced = 1; %make sure advanced setting is on
                    disp_identification(pdraws(jcrit(j),:), ide_reducedform_(j), ide_moments_(j), ide_spectrum_(j), ide_minimal_(j), name, options_ident);
                    options_ident.advanced = advanced0; % reset advanced                    
                    if ~options_.nograph
                        % plot (i) identification strength and sensitivity measure based on the sample information matrix and (ii) advanced analysis graphs
                        plot_identification(pdraws(jcrit(j),:), ide_moments_(j), ide_hess_(j), ide_reducedform_(j), ide_lre_(j), 1, tittxt, name, IdentifDirectoryName, tittxt, name_tex);
                    end
                end
                if ~iload
                    save([IdentifDirectoryName '/' fname '_identif.mat'], 'ide_hess_', 'ide_moments_', 'ide_reducedform_', 'ide_lre_', 'ide_spectrum_', 'ide_minimal_', 'jcrit', '-append');
                end
                % reset nodisplay option
                options_.nodisplay = store_nodisplay;                
            end
        end
    end
end

%reset warning state
if isoctave
    warning('on')
else
    warning on
end

skipline()
disp('==== Identification analysis completed ====')
skipline(2)

options_ = store_options_; %restore options set
