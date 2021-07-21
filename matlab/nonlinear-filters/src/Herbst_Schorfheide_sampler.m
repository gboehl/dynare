function Herbst_Schorfheide_sampler(TargetFun,xparam1,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% function Herbst_Schorfheide_sampler(TargetFun,xparam1,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% SMC sampler from JAE 2014 .
%
% INPUTS
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters.
%   o dataset_              data structure
%   o dataset_info          dataset info structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code suitable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle)
% has been removed from this function and been placed in the posterior_sampler_core.m funtion.
%
% The DYNARE parallel packages comprise a i) set of pairs of Matlab functions that can be executed in
% parallel and called name_function.m and name_function_core.m and ii) a second set of functions used
% to manage the parallel computations.
%
% This function was the first function to be parallelized. Later, other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management functions.

% Copyright (C) 2006-2021 Dynare Team
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


% Create the tempering schedule
phi = bsxfun(@power,(bsxfun(@minus,1:1:options_.HSsmc.nphi,1)/(options_.HSsmc.nphi-1)),options_.HSsmc.lambda) ;
% tuning for MH algorithms matrices
zhat    = 0 ;                                    % normalization constant
csim    = zeros(options_.HSsmc.nphi,1) ;         % scale parameter
ESSsim  = zeros(options_.HSsmc.nphi,1) ;         % ESS
acptsim = zeros(options_.HSsmc.nphi,1) ;         % average acceptance rate
% Step 0: Initialization of the sampler
[ param, tlogpost_i, loglik, npar, ~, bayestopt_] = ...
    SMC_samplers_initialization(TargetFun, xparam1, mh_bounds, dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,options_.HSsmc.nparticles);
weights = ones(options_.HSsmc.nparticles,1)/options_.HSsmc.nparticles ;
% The Herbst and Schorfheide sampler starts here
for i=2:options_.HSsmc.nphi
    % (a) Correction
    % incremental weights
    incwt = exp((phi(i)-phi(i-1))*loglik) ;
    % update weights
    weights = bsxfun(@times,weights,incwt) ;
    sum_weights = sum(weights) ;
    zhat = zhat + log(sum_weights) ;
    % normalize weights
    weights = weights/sum_weights ;
    % (b) Selection
    ESSsim(i) = 1/sum(weights.^2) ;
    if (ESSsim(i) < options_.HSsmc.nparticles/2)
        indx_resmpl = smc_resampling(weights,rand(1,1),options_.HSsmc.nparticles) ;
        param = param(:,indx_resmpl) ;
        loglik = loglik(indx_resmpl) ;
        tlogpost_i = tlogpost_i(indx_resmpl) ;
        weights = ones(options_.HSsmc.nparticles,1)/options_.HSsmc.nparticles  ;
    end
    % (c) Mutation
    options_.HSsmc.c = options_.HSsmc.c*modified_logit(0.95,0.1,16.0,options_.HSsmc.acpt-options_.HSsmc.trgt) ;
    % Calculate estimates of mean and variance
    mu = param*weights ;
    z = bsxfun(@minus,param,mu) ;
    R = z*(bsxfun(@times,z',weights)) ;
    Rchol = chol(R)' ;
    % Mutation
    if options_.HSsmc.option_mutation==1
        [param,tlogpost_i,loglik,options_.HSsmc.acpt] = mutation_RW(TargetFun,param,tlogpost_i,loglik,phi,i,options_.HSsmc.c*Rchol,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_) ;
    elseif options_.HSsmc.option_mutation==2
        inv_R = inv(options_.HSsmc.c^2*R) ;
        Rdiagchol = sqrt(diag(R)) ;
        [param,tlogpost_i,loglik,options_.HSsmc.acpt] = mutation_Mixture(TargetFun,param,tlogpost_i,loglik,phi,i,options_.HSsmc.c*Rchol,options_.HSsmc.c*Rdiagchol,inv_R,mu,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_) ;
    end
    acptsim(i) = options_.HSsmc.acpt ;
    csim(i) = options_.HSsmc.c ;
    % print information
    fprintf(' Iteration = %5.0f / %5.0f \n', i, options_.HSsmc.nphi);
    fprintf(' phi = %5.4f \n', phi(i));
    fprintf(' Neff = %5.4f \n', ESSsim(i));
    fprintf(' %accept. = %5.4f \n', acptsim(i));
end
indx_resmpl = smc_resampling(weights,rand(1,1),options_.HSsmc.nparticles);
distrib_param = param(:,indx_resmpl);
fprintf(' Log_lik = %5.4f \n', zhat);

mean_xparam = mean(distrib_param,2);
%mat_var_cov = bsxfun(@minus,distrib_param,mean_xparam) ;
%mat_var_cov = (mat_var_cov*mat_var_cov')/(options_.HSsmc.nparticles-1) ;
%std_xparam = sqrt(diag(mat_var_cov)) ;
lb95_xparam = zeros(npar,1) ;
ub95_xparam = zeros(npar,1) ;
for i=1:npar
    temp = sortrows(distrib_param(i,:)') ;
    lb95_xparam(i) = temp(0.025*options_.HSsmc.nparticles) ;
    ub95_xparam(i) = temp(0.975*options_.HSsmc.nparticles) ;
end

TeX = options_.TeX;

str = sprintf(' Param. \t Lower Bound (95%%) \t Mean \t Upper Bound (95%%)');
for l=1:npar
    [name,~] = get_the_name(l,TeX,M_,estim_params_,options_);
    str = sprintf('%s\n %s \t\t %5.4f \t\t %7.5f \t\t %5.4f', str, name, lb95_xparam(l), mean_xparam(l), ub95_xparam(l));
end
disp([str])
disp('')

%% Plot parameters densities

[nbplt,nr,nc,lr,lc,nstar] = pltorg(npar);

if TeX
    fidTeX = fopen([M_.fname '_param_density.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by DSMH.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourier Transform approximation.
plt = 1 ;
%for plt = 1:nbplt,
if TeX
    NAMES = [];
    TeXNAMES = [];
end
hh = dyn_figure(options_.nodisplay,'Name','Parameters Densities');
for k=1:npar %min(nstar,npar-(plt-1)*nstar)
    subplot(ceil(sqrt(npar)),floor(sqrt(npar)),k)
    %kk = (plt-1)*nstar+k;
    [name,texname] = get_the_name(k,TeX,M_,estim_params_,options_);
    optimal_bandwidth = mh_optimal_bandwidth(distrib_param(k,:)',options_.HSsmc.nparticles,bandwidth,kernel_function);
    [density(:,1),density(:,2)] = kernel_density_estimate(distrib_param(k,:)',number_of_grid_points,...
        options_.HSsmc.nparticles,optimal_bandwidth,kernel_function);
    plot(density(:,1),density(:,2));
    hold on
    if TeX
        title(texname,'interpreter','latex')
    else
        title(name,'interpreter','none')
    end
    hold off
    axis tight
    drawnow
end
dyn_saveas(hh,[ M_.fname '_param_density' int2str(plt) ],options_.nodisplay,options_.graph_format);
if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    % TeX eps loader file
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%_param_density%s}\n',min(k/floor(sqrt(npar)),1),M_.fname,int2str(plt));
    fprintf(fidTeX,'\\caption{Parameter densities based on the Herbst/Schorfheide sampler.}');
    fprintf(fidTeX,'\\label{Fig:ParametersDensities:%s}\n',int2str(plt));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end
%end

function [param,tlogpost_i,loglik,acpt] = mutation_RW(TargetFun,param,tlogpost_i,loglik,phi,i,cRchol,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
acpt = 0.0 ;
tlogpost_i = tlogpost_i + (phi(i)-phi(i-1))*loglik ;
for j=1:options_.HSsmc.nparticles
    validate= 0;
    while validate==0
        candidate = param(:,j) + cRchol*randn(size(param,1),1) ;
        if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
            [tlogpostx,loglikx] = tempered_likelihood(TargetFun,candidate,phi(i),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
            if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                validate = 1;
                if rand(1,1)<exp(tlogpostx-tlogpost_i(j)) % accept
                    acpt = acpt + 1 ;
                    param(:,j) = candidate;
                    loglik(j) = loglikx;
                    tlogpost_i(j) = tlogpostx;
                end
            end
        end
    end
end
acpt = acpt/options_.HSsmc.nparticles;

function [param,tlogpost_i,loglik,acpt] = mutation_Mixture(TargetFun,param,tlogpost_i,loglik,phi,i,cRchol,cRdiagchol,invR,mu,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
acpt = 0.0 ;
tlogpost_i = tlogpost_i + (phi(i)-phi(i-1))*loglik ;
for j=1:options_.HSsmc.nparticles
    qx = 0 ;
    q0 = 0 ;
    alpind = rand(1,1) ;
    validate= 0;
    while validate==0
        if alpind<options_.HSsmc.alp % RW, no need to modify qx and q0
            candidate = param(:,j) + cRchol*randn(size(param,1),1);
        elseif alpind<options_.HSsmc.alp + (1-options_.HSsmc.alp/2) % random walk with diagonal cov no need to modify qx and q0
            candidate = param(:,j) + cRdiagchol*randn(size(param,1),1);
        else % Proposal densities
            candidate = mu + cRchol*randn(size(param,1),1);
            qx = -.5*(candidate-mu)'*invR*(candidate-mu) ; % no need of the constants in the acceptation rule
            q0 = -.5*(param(:,j)-mu)'*invR*(param(:,j)-mu) ;
        end
        if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
            [tlogpostx,loglikx] = tempered_likelihood(TargetFun,candidate,phi(i),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
            if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                validate = 1;
                if rand(1,1)<exp((tlogpostx-qx)-(tlogpost_i(j)-q0)) % accept
                    acpt = acpt + 1 ;
                    param(:,j) = candidate;
                    loglik(j) = loglikx;
                    tlogpost_i(j) = tlogpostx;
                end
            end
        end
    end
end
acpt = acpt/options_.HSsmc.nparticles;

function x = modified_logit(a,b,c,d)
tmp = exp(c*d) ;
x =  a + b*tmp/(1 + tmp) ;
