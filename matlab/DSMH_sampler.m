function DSMH_sampler(TargetFun,xparam1,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% function DSMH_sampler(TargetFun,xparam1,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% Dynamic Striated Metropolis-Hastings algorithm.
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

% Copyright © 2006-2022 Dynare Team
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


lambda = exp(bsxfun(@minus,options_.posterior_sampler_options.dsmh.H,1:1:options_.posterior_sampler_options.dsmh.H)/(options_.posterior_sampler_options.dsmh.H-1)*log(options_.posterior_sampler_options.dsmh.lambda1));
c = 0.055 ;
MM = int64(options_.posterior_sampler_options.dsmh.N*options_.posterior_sampler_options.dsmh.G/10) ;

% Step 0: Initialization of the sampler
[ param, tlogpost_iminus1, loglik, bayestopt_] = ...
    SMC_samplers_initialization(TargetFun, xparam1, mh_bounds, dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,options_.posterior_sampler_options.dsmh.nparticles);

ESS = zeros(options_.posterior_sampler_options.dsmh.H,1) ;
zhat = 1 ;

% The DSMH starts here
for i=2:options_.posterior_sampler_options.dsmh.H
    disp('');
    disp('Tempered iteration');
    disp(i) ;
    % Step 1: sort the densities and compute IS weigths
    [tlogpost_iminus1,loglik,param] = sort_matrices(tlogpost_iminus1,loglik,param) ;
    [tlogpost_i,weights,zhat,ESS,Omegachol] = compute_IS_weights_and_moments(param,tlogpost_iminus1,loglik,lambda,i,zhat,ESS) ;
    % Step 2: tune c_i
    c = tune_c(TargetFun,param,tlogpost_i,lambda,i,c,Omegachol,weights,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
    % Step 3: Metropolis step
    [param,tlogpost_iminus1,loglik] = mutation_DSMH(TargetFun,param,tlogpost_i,tlogpost_iminus1,loglik,lambda,i,c,MM,Omegachol,weights,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
end

weights = exp(loglik*(lambda(end)-lambda(end-1)));
weights = weights/sum(weights);
indx_resmpl = smc_resampling(weights,rand(1,1),options_.posterior_sampler_options.dsmh.nparticles);
distrib_param = param(:,indx_resmpl);

mean_xparam = mean(distrib_param,2);
npar  = length(xparam1);
%mat_var_cov = bsxfun(@minus,distrib_param,mean_xparam) ;
%mat_var_cov = (mat_var_cov*mat_var_cov')/(options_.HSsmc.nparticles-1) ;
%std_xparam = sqrt(diag(mat_var_cov)) ;
lb95_xparam = zeros(npar,1) ;
ub95_xparam = zeros(npar,1) ;
for i=1:npar
    temp = sortrows(distrib_param(i,:)') ;
    lb95_xparam(i) = temp(0.025*options_.posterior_sampler_options.dsmh.nparticles) ;
    ub95_xparam(i) = temp(0.975*options_.posterior_sampler_options.dsmh.nparticles) ;
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
    optimal_bandwidth = mh_optimal_bandwidth(distrib_param(k,:)',options_.posterior_sampler_options.dsmh.nparticles,bandwidth,kernel_function);
    [density(:,1),density(:,2)] = kernel_density_estimate(distrib_param(k,:)',number_of_grid_points,...
        options_.posterior_sampler_options.dsmh.nparticles,optimal_bandwidth,kernel_function);
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
    fprintf(fidTeX,'\\caption{Parameter densities based on the Dynamic Striated Metropolis-Hastings algorithm.}');
    fprintf(fidTeX,'\\label{Fig:ParametersDensities:%s}\n',int2str(plt));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end
%end

function [tlogpost_iminus1,loglik,param] = sort_matrices(tlogpost_iminus1,loglik,param)
[~,indx_ord] = sortrows(tlogpost_iminus1);
tlogpost_iminus1 = tlogpost_iminus1(indx_ord);
param = param(:,indx_ord);
loglik = loglik(indx_ord);

function [tlogpost_i,weights,zhat,ESS,Omegachol] = compute_IS_weights_and_moments(param,tlogpost_iminus1,loglik,lambda,i,zhat,ESS)
if i==1
    tlogpost_i = tlogpost_iminus1 + loglik*lambda(i);
else
    tlogpost_i = tlogpost_iminus1 + loglik*(lambda(i)-lambda(i-1));
end
weights = exp(tlogpost_i-tlogpost_iminus1);
zhat = (mean(weights))*zhat ;
weights = weights/sum(weights);
ESS(i) = 1/sum(weights.^2);
% estimates of mean and variance
mu = param*weights;
z = bsxfun(@minus,param,mu);
Omega = z*diag(weights)*z';
Omegachol = chol(Omega)';

function c = tune_c(TargetFun,param,tlogpost_i,lambda,i,c,Omegachol,weights,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
disp('tuning c_i...');
disp('Initial value =');
disp(c) ;
npar = size(param,1);
lower_prob = (.5*(options_.posterior_sampler_options.dsmh.alpha0+options_.posterior_sampler_options.dsmh.alpha1))^5;
upper_prob = (.5*(options_.posterior_sampler_options.dsmh.alpha0+options_.posterior_sampler_options.dsmh.alpha1))^(1/5);
stop=0 ;
while stop==0
    acpt = 0.0;
    indx_resmpl = smc_resampling(weights,rand(1,1),options_.posterior_sampler_options.dsmh.G);
    param0 = param(:,indx_resmpl);
    tlogpost0 = tlogpost_i(indx_resmpl);
    for j=1:options_.posterior_sampler_options.dsmh.G
        for l=1:options_.posterior_sampler_options.dsmh.K
            validate = 0;
            while validate == 0
                candidate = param0(:,j) + sqrt(c)*Omegachol*randn(npar,1);
                if all(candidate >= mh_bounds.lb) && all(candidate <= mh_bounds.ub)
                    [tlogpostx,loglikx] = tempered_likelihood(TargetFun,candidate,lambda(i),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
                    if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                        validate = 1;
                        if rand(1,1)<exp(tlogpostx-tlogpost0(j)) % accept
                            acpt = acpt + 1/(options_.posterior_sampler_options.dsmh.G*options_.posterior_sampler_options.dsmh.K);
                            param0(:,j)= candidate;
                            tlogpost0(j) = tlogpostx;
                        end
                    end
                end
            end
        end
    end
    disp('Acceptation rate =') ;
    disp(acpt) ;
    if options_.posterior_sampler_options.dsmh.alpha0<=acpt && acpt<=options_.posterior_sampler_options.dsmh.alpha1
        disp('done!');
        stop=1;
    else
        if acpt<lower_prob
            c = c/5;
        elseif lower_prob<=acpt && acpt<=upper_prob
            c = c*log(.5*(options_.posterior_sampler_options.dsmh.alpha0+options_.posterior_sampler_options.dsmh.alpha1))/log(acpt);
        else
            c = 5*c;
        end
        disp('Trying with c= ') ;
        disp(c)
    end
end

function [out_param,out_tlogpost_iminus1,out_loglik] = mutation_DSMH(TargetFun,param,tlogpost_i,tlogpost_iminus1,loglik,lambda,i,c,MM,Omegachol,weights,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
indx_levels = (1:1:MM-1)*options_.posterior_sampler_options.dsmh.N*options_.posterior_sampler_options.dsmh.G/MM;
npar = size(param,1) ;
p = 1/(10*options_.posterior_sampler_options.dsmh.tau);
disp('Metropolis step...');
% build the dynamic grid of levels
levels = [0.0;tlogpost_iminus1(indx_levels)];
% initialize the outputs
out_param = param;
out_tlogpost_iminus1 = tlogpost_i;
out_loglik = loglik;
% resample and initialize the starting groups
indx_resmpl = smc_resampling(weights,rand(1,1),options_.posterior_sampler_options.dsmh.G);
param0 = param(:,indx_resmpl);
tlogpost_iminus10 = tlogpost_iminus1(indx_resmpl);
tlogpost_i0 = tlogpost_i(indx_resmpl);
loglik0 = loglik(indx_resmpl);
% Start the Metropolis
for l=1:options_.posterior_sampler_options.dsmh.N*options_.posterior_sampler_options.dsmh.tau
    for j=1:options_.posterior_sampler_options.dsmh.G
        u1 = rand(1,1);
        u2 = rand(1,1);
        if u1<p
            k=1 ;
            for m=1:MM-1
                if levels(m)<=tlogpost_iminus10(j) && tlogpost_iminus10(j)<levels(m+1)
                    k = m+1;
                    break
                end
            end
            indx = floor( (k-1)*options_.posterior_sampler_options.dsmh.N*options_.posterior_sampler_options.dsmh.G/MM+1 + u2*(options_.posterior_sampler_options.dsmh.N*options_.posterior_sampler_options.dsmh.G/MM-1) );
            if i==1
                alp = (loglik(indx)-loglik0(j))*lambda(i);
            else
                alp = (loglik(indx)-loglik0(j))*(lambda(i)-lambda(i-1));
            end
            if u2<exp(alp)
                param0(:,j) = param(:,indx);
                tlogpost_i0(j) = tlogpost_i(indx);
                loglik0(j) = loglik(indx);
                tlogpost_iminus10(j) = tlogpost_iminus1(indx);
            end
        else
            validate= 0;
            while validate==0
                candidate = param0(:,j) + sqrt(c)*Omegachol*randn(npar,1);
                if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
                    [tlogpostx,loglikx] = tempered_likelihood(TargetFun,candidate,lambda(i),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
                    if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                        validate = 1;
                        if u2<exp(tlogpostx-tlogpost_i0(j)) % accept
                            param0(:,j) = candidate;
                            tlogpost_i0(j) = tlogpostx;
                            loglik0(j) = loglikx;
                            if i==1
                                tlogpost_iminus10(j) = tlogpostx-loglikx*lambda(i);
                            else
                                tlogpost_iminus10(j) = tlogpostx-loglikx*(lambda(i)-lambda(i-1));
                            end
                        end
                    end
                end
            end
        end
    end
    if mod(l,options_.posterior_sampler_options.dsmh.tau)==0
        rang = (l/options_.posterior_sampler_options.dsmh.tau-1)*options_.posterior_sampler_options.dsmh.G+1:l*options_.posterior_sampler_options.dsmh.G/options_.posterior_sampler_options.dsmh.tau;
        out_param(:,rang) = param0;
        out_tlogpost_iminus1(rang) = tlogpost_i0;
        out_loglik(rang) = loglik0;
    end
end