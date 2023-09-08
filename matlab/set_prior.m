function [xparam1, estim_params_, bayestopt_, lb, ub, M_]=set_prior(estim_params_, M_, options_)
% function [xparam1,estim_params_,bayestopt_,lb,ub, M_]=set_prior(estim_params_, M_, options_)
% sets prior distributions
%
% INPUTS
%    o estim_params_    [structure] characterizing parameters to be estimated.
%    o M_               [structure] characterizing the model.
%    o options_         [structure] characterizing the options.
%
% OUTPUTS
%    o xparam1          [double]    vector of parameters to be estimated (initial values)
%    o estim_params_    [structure] characterizing parameters to be estimated
%    o bayestopt_       [structure] characterizing priors
%    o lb               [double]    vector of lower bounds for the estimated parameters.
%    o ub               [double]    vector of upper bounds for the estimated parameters.
%    o M_               [structure] characterizing the model.
%
% SPECIAL REQUIREMENTS
%    None

% Copyright © 2003-2023 Dynare Team
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

nvx = size(estim_params_.var_exo,1);
nvn = size(estim_params_.var_endo,1);
ncx = size(estim_params_.corrx,1);
ncn = size(estim_params_.corrn,1);
np = size(estim_params_.param_vals,1);

estim_params_.nvx = nvx; %exogenous shock variances
estim_params_.nvn = nvn; %endogenous variances, i.e. measurement error
estim_params_.ncx = ncx; %exogenous shock correlations
estim_params_.ncn = ncn; % correlation between endogenous variables, i.e. measurement error.
estim_params_.np = np;   % other parameters of the model

xparam1 = [];
ub = []; % Upper bound imposed for optimization.
lb = []; % Lower bound imposed for optimization.
bayestopt_.pshape = [];
bayestopt_.p1 = []; % prior mean
bayestopt_.p2 = []; % prior standard deviation
bayestopt_.p3 = []; % lower bound of the distribution, only considering whether a generalized distribution is used, not when the prior is truncated
bayestopt_.p4 = []; % upper bound of the distribution, only considering whether a generalized distribution is used, not when the prior is truncated
bayestopt_.p5 = zeros(nvx+nvn+ncx+ncn+np,1); % prior mode
bayestopt_.p6 = []; % first hyper-parameter (\alpha for the BETA and GAMMA distributions, s for the INVERSE GAMMAs, expectation for the GAUSSIAN distribution, lower bound for the UNIFORM distribution).
bayestopt_.p7 = []; % second hyper-parameter (\beta for the BETA and GAMMA distributions, \nu for the INVERSE GAMMAs, standard deviation for the GAUSSIAN distribution, upper bound for the UNIFORM distribution).

bayestopt_.jscale = []; %jscale might subsequently be overwritten by mode_compute=6 or check_posterior_sampler_options
bayestopt_.name = {};
if nvx
    xparam1 = estim_params_.var_exo(:,2);
    ub = estim_params_.var_exo(:,4);
    lb = estim_params_.var_exo(:,3);
    bayestopt_.pshape =  estim_params_.var_exo(:,5);
    bayestopt_.p1 =  estim_params_.var_exo(:,6);
    bayestopt_.p2 =  estim_params_.var_exo(:,7);
    bayestopt_.p3 =  estim_params_.var_exo(:,8); %take generalized distribution into account
    bayestopt_.p4 =  estim_params_.var_exo(:,9); %take generalized distribution into account
    bayestopt_.jscale =  estim_params_.var_exo(:,10);
    bayestopt_.name = M_.exo_names(estim_params_.var_exo(:,1));
end
if nvn
    estim_params_.nvn_observable_correspondence=NaN(nvn,1); % stores number of corresponding observable
    if isequal(M_.H,0) %if no previously set measurement error, initialize H
        nvarobs = length(options_.varobs);
        M_.H = zeros(nvarobs,nvarobs);
        M_.Correlation_matrix_ME = eye(nvarobs);
    end
    for i=1:nvn
        obsi_ = strmatch(M_.endo_names{estim_params_.var_endo(i,1)}, options_.varobs, 'exact');
        if isempty(obsi_)
            error(['The variable ' M_.endo_names{estim_params_.var_endo(i,1)} ' has to be declared as observable since you assume a measurement error on it.'])
        end
        estim_params_.nvn_observable_correspondence(i,1)=obsi_;
    end
    xparam1 = [xparam1; estim_params_.var_endo(:,2)];
    ub = [ub; estim_params_.var_endo(:,4)];
    lb = [lb; estim_params_.var_endo(:,3)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.var_endo(:,5)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.var_endo(:,6)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.var_endo(:,7)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.var_endo(:,8)]; %take generalized distribution into account
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.var_endo(:,9)]; %take generalized distribution into account
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.var_endo(:,10)];
    bayestopt_.name = [ bayestopt_.name; options_.varobs(estim_params_.nvn_observable_correspondence)];
end
if ncx
    xparam1 = [xparam1; estim_params_.corrx(:,3)];
    ub = [ub; max(min(estim_params_.corrx(:,5),1),-1)];
    lb = [lb; min(max(estim_params_.corrx(:,4),-1),1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrx(:,6)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.corrx(:,7)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.corrx(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrx(:,9)]; %take generalized distribution into account
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrx(:,10)]; %take generalized distribution into account
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrx(:,11)];
    baseid = length(bayestopt_.name);
    bayestopt_.name = [bayestopt_.name; cell(ncx, 1)];
    for i = 1:ncx
        bayestopt_.name(baseid+i) = {sprintf('corr %s, %s', ...
                                             M_.exo_names{estim_params_.corrx(i,1)}, ...
                                             M_.exo_names{estim_params_.corrx(i,2)})};
    end
end
if ncn
    estim_params_.corrn_observable_correspondence=NaN(ncn,2);
    if isequal(M_.H,0)
        nvarobs = length(options_.varobs);
        M_.H = zeros(nvarobs,nvarobs);
        M_.Correlation_matrix_ME = eye(nvarobs);
    end
    xparam1 = [xparam1; estim_params_.corrn(:,3)];
    ub = [ub; max(min(estim_params_.corrn(:,5),1),-1)];
    lb = [lb; min(max(estim_params_.corrn(:,4),-1),1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrn(:,6)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.corrn(:,7)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.corrn(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrn(:,9)]; %take generalized distribution into account
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrn(:,10)]; %take generalized distribution into account
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrn(:,11)];
    baseid = length(bayestopt_.name);
    bayestopt_.name = [bayestopt_.name; cell(ncn, 1)];;
    for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        bayestopt_.name(baseid+i) = {sprintf('corr %s, %s', M_.endo_names{k1}, M_.endo_names{k2})};
        % find correspondence to varobs to construct H in set_all_paramters
        obsi1 = strmatch(M_.endo_names{k1}, options_.varobs, 'exact');
        obsi2 = strmatch(M_.endo_names{k2}, options_.varobs, 'exact');
        % save correspondence
        estim_params_.corrn_observable_correspondence(i,:)=[obsi1, obsi2];
    end
end
if np
    xparam1 = [xparam1; estim_params_.param_vals(:,2)];
    ub = [ub; estim_params_.param_vals(:,4)];
    lb = [lb; estim_params_.param_vals(:,3)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.param_vals(:,5)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.param_vals(:,6)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.param_vals(:,7)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.param_vals(:,8)]; %take generalized distribution into account
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.param_vals(:,9)]; %take generalized distribution into account
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.param_vals(:,10)];
    bayestopt_.name = [bayestopt_.name; M_.param_names(estim_params_.param_vals(:,1))];
end

bayestopt_.p6 = NaN(size(bayestopt_.p1)) ;
bayestopt_.p7 = bayestopt_.p6 ;

%% check for point priors and disallow them as they do not work with MCMC
if any(bayestopt_.p2 ==0)
    error(sprintf(['Error in prior for %s: you cannot use a point prior in estimation. Either increase the prior standard deviation',...
                   ' or fix the parameter completely.'], bayestopt_.name{bayestopt_.p2 ==0}))
end

% generalized location parameters by default for beta distribution
k = find(bayestopt_.pshape == 1);
k1 = find(isnan(bayestopt_.p3(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
k1 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p4(k(k1)) = ones(length(k1),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)), bayestopt_.p7(k(i))] = beta_specification(bayestopt_.p1(k(i)), bayestopt_.p2(k(i))^2, bayestopt_.p3(k(i)), bayestopt_.p4(k(i)), bayestopt_.name{k(i)});
    if bayestopt_.p6(k(i))<1 || bayestopt_.p7(k(i))<1
        fprintf('Prior distribution for parameter %s has unbounded density!\n',bayestopt_.name{k(i)})
    end
    m = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) , bayestopt_.p4(k(i)) ],1);
    if length(m)==1
        bayestopt_.p5(k(i)) = m;
    else
        disp(['Prior distribution for parameter ' bayestopt_.name{k(i)}  ' has two modes!'])
        bayestopt_.p5(k(i)) = m(1);
    end
end

% generalized location parameter by default for gamma distribution
k =  find(bayestopt_.pshape == 2);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)), bayestopt_.p7(k(i))] = gamma_specification(bayestopt_.p1(k(i)), bayestopt_.p2(k(i))^2, bayestopt_.p3(k(i)), bayestopt_.name{k(i)});
    if bayestopt_.p6(k(i))<1
        fprintf('Prior distribution for parameter %s has unbounded density!\n',bayestopt_.name{k(i)})
    end
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 2) ;
end

% truncation parameters by default for normal distribution
k  = find(bayestopt_.pshape == 3);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = -Inf*ones(length(k1),1);
bayestopt_.p4(k(k2)) = Inf*ones(length(k2),1);
bayestopt_.p6(k) = bayestopt_.p1(k);
bayestopt_.p7(k) = bayestopt_.p2(k);
bayestopt_.p5(k) = bayestopt_.p1(k);

% inverse gamma distribution (type 1)
k = find(bayestopt_.pshape == 4);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = inverse_gamma_specification(bayestopt_.p1(k(i)), bayestopt_.p2(k(i))^2, bayestopt_.p3(k(i)), 1, false, bayestopt_.name{k(i)});
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 4);
end

% uniform distribution
k = find(bayestopt_.pshape == 5);
problem_parameters='';
for i=1:length(k)
    [bayestopt_.p1(k(i)),bayestopt_.p2(k(i)),bayestopt_.p6(k(i)),bayestopt_.p7(k(i)),error_indicator] = ...
        uniform_specification(bayestopt_.p1(k(i)),bayestopt_.p2(k(i)),bayestopt_.p3(k(i)),bayestopt_.p4(k(i)));
    if error_indicator
        if isempty(problem_parameters)
            problem_parameters=[bayestopt_.name{k(i)}];
        else
            problem_parameters=[problem_parameters ', ' bayestopt_.name{k(i)}];
        end
    end
    bayestopt_.p3(k(i)) = bayestopt_.p6(k(i)) ;
    bayestopt_.p4(k(i)) = bayestopt_.p7(k(i)) ;
    bayestopt_.p5(k(i)) = NaN ;
end
if ~isempty(problem_parameters)
    error(['uniform_specification: You defined lower and upper bounds for parameter ', problem_parameters, '. In this case, you need to leave mean and standard deviation empty.'])
end

% inverse gamma distribution (type 2)
k = find(bayestopt_.pshape == 6);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = ...
        inverse_gamma_specification(bayestopt_.p1(k(i)), bayestopt_.p2(k(i))^2, bayestopt_.p3(k(i)), 2, false, bayestopt_.name{k(i)});
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 6) ;
end

% Weibull distribution
k = find(bayestopt_.pshape == 8);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = weibull_specification(bayestopt_.p1(k(i)), bayestopt_.p2(k(i))^2, bayestopt_.p3(k(i)), bayestopt_.name{k(i)});
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 8) ;
end

k = find(isnan(xparam1));
if ~isempty(k)
    xparam1(k) = bayestopt_.p1(k);
end

if options_.initialize_estimated_parameters_with_the_prior_mode
    xparam1 = bayestopt_.p5;
    k = find(isnan(xparam1));% Because the uniform density do not have a mode!
    if ~isempty(k)
        xparam1(k) = bayestopt_.p1(k);
    end
end



% I create subfolder M_.dname/prior if needed.
CheckPath('prior',M_.dname);

% I save the prior definition if the prior has changed.
if exist([ M_.dname '/prior/definition.mat'])
    old = load([M_.dname '/prior/definition.mat'],'bayestopt_');
    prior_has_changed = 0;
    if length(bayestopt_.p1)==length(old.bayestopt_.p1)
        if any(bayestopt_.p1-old.bayestopt_.p1)
            prior_has_changed = 1;
        elseif any(bayestopt_.p2-old.bayestopt_.p2)
            prior_has_changed = 1;
        elseif any(bayestopt_.p3-old.bayestopt_.p3)
            prior_has_changed = 1;
        elseif any(bayestopt_.p4-old.bayestopt_.p4)
            prior_has_changed = 1;
        elseif any(bayestopt_.p5-old.bayestopt_.p5(:))
            prior_has_changed = 1;
        elseif any(bayestopt_.p6-old.bayestopt_.p6)
            prior_has_changed = 1;
        elseif any(bayestopt_.p7-old.bayestopt_.p7)
            prior_has_changed = 1;
        end
    else
        prior_has_changed = 1;
    end
    if prior_has_changed
        delete([M_.dname '/prior/definition.mat']);
        save([M_.dname '/prior/definition.mat'],'bayestopt_');
    end
else
    save([M_.dname '/prior/definition.mat'],'bayestopt_');
end

% initialize persistent variables in priordens()
priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, ...
          bayestopt_.p3,bayestopt_.p4,1);