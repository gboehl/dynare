function [ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags)
%function [ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags)
% bvar_toolbox  Routines shared between BVAR methods
% Computes several things for the estimations of a BVAR(nlags)
%
% INPUTS:
%     nlags: number of lags
%
% OUTPUTS:
%    ny:            number of endogenous variables
%    nx:            number of exogenous variables (equal to zero, or one if a
%                   constant term is included)
%    posterior:     a structure describing the posterior distribution (which is
%                   normal-Inverse-Wishart)
%                   Its fields are:
%                   - df: degrees of freedom of the inverse-Wishart distribution
%                   - S: matrix parameter for the inverse-Wishart distribution
%                   - XXi: first component of the VCV of the matrix-normal
%                     distribution (the other one being drawn from the
%                     inverse-Wishart)
%                   - PhiHat: mean of the matrix-normal distribution
%    prior:         a structure describing the prior distribution
%                   Its fields are the same than for the posterior
%    forecast_data: a structure containing data useful for forecasting
%                   Its fields are:
%                   - initval: a nlags*ny matrix containing the "nlags" last
%                     observations of the sample (i.e. before options_.nobs)
%                   - xdata: a matrix containing the future exogenous for
%                     forecasting, of size options_.forecast*nx (actually only
%                     contains "1" values for the constant term if nx ~= 0)
%                   - realized_val: only non-empty if options_.nobs doesn't point
%                     to the end of sample    
%                     In that case, contains values of endogenous variables after
%                     options_.nobs and up to the end of the sample
%                   - realized_xdata: contains values of exogenous variables after
%                     options_.nobs and up to the end of the sample (actually only
%                     contains "1" values for the constant term if nx ~= 0)
%
% SPECIAL REQUIREMENTS:
%    This function uses the following Dynare options:
%    - datafile, first_obs, varobs, xls_sheet, xls_range, nobs, presample
%    - bvar_prior_{tau,decay,lambda,mu,omega,flat,train}

% Copyright (C) 2003-2007 Christopher Sims
% Copyright (C) 2007-2012 Dynare Team
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

global options_

% Load dataset
dataset = read_variables(options_.datafile, options_.varobs, [], options_.xls_sheet, options_.xls_range);
options_ = set_default_option(options_, 'nobs', size(dataset,1)-options_.first_obs+1);

if (options_.first_obs+options_.nobs-1)> size(dataset,1)
    fprintf('Incorrect or missing specification of the number of observations. nobs can be at most %4u\n',size(dataset,1)-options_.first_obs+1);
    error('Inconsistent number of observations.') 
end

% Parameters for prior
if options_.first_obs + options_.presample <= nlags
    error('first_obs+presample should be > nlags (for initializing the VAR)')
end

train = options_.bvar_prior_train;

if options_.first_obs + options_.presample - train <= nlags
    error('first_obs+presample-train should be > nlags (for initializating the VAR)')
end

idx = options_.first_obs+options_.presample-train-nlags:options_.first_obs+options_.nobs-1;

% Prepare dataset
if options_.loglinear && ~options_.logdata
    dataset = log(dataset);
end
if options_.prefilter
    dataset(idx,:) = dataset(idx,:) - ones(length(idx),1)*mean(dataset(idx,:));
end

mnprior.tight = options_.bvar_prior_tau;
mnprior.decay = options_.bvar_prior_decay;

% Use only initializations lags for the variance prior
vprior.sig = std(dataset(options_.first_obs+options_.presample-nlags:options_.first_obs+options_.presample,:))';
vprior.w = options_.bvar_prior_omega;

lambda = options_.bvar_prior_lambda;
mu = options_.bvar_prior_mu;
flat = options_.bvar_prior_flat;

ny = size(dataset, 2);
if options_.prefilter || options_.noconstant
    nx = 0;
else
    nx = 1;
end

[ydum, xdum, pbreaks] = varprior(ny, nx, nlags, mnprior, vprior);

ydata = dataset(idx, :);
T = size(ydata, 1);
xdata = ones(T,nx);

% Posterior density
var = rfvar3([ydata; ydum], nlags, [xdata; xdum], [T; T+pbreaks], lambda, mu);
Tu = size(var.u, 1);

posterior.df = Tu - ny*nlags - nx - flat*(ny+1);
posterior.S = var.u' * var.u;
posterior.XXi = var.xxi;
posterior.PhiHat = var.B;

% Prior density
Tp = train + nlags;
if nx
    xdata = xdata(1:Tp, :);
else
    xdata = [];
end
varp = rfvar3([ydata(1:Tp, :); ydum], nlags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
Tup = size(varp.u, 1);

prior.df = Tup - ny*nlags - nx - flat*(ny+1);
prior.S = varp.u' * varp.u;
prior.XXi = varp.xxi;
prior.PhiHat = varp.B;

if prior.df < ny
    error('Too few degrees of freedom in the inverse-Wishart part of prior distribution. You should increase training sample size.')
end

% Add forecast informations
if nargout >= 5
    forecast_data.xdata = ones(options_.forecast, nx);
    forecast_data.initval = ydata(end-nlags+1:end, :);
    if options_.first_obs + options_.nobs <= size(dataset, 1)
        forecast_data.realized_val = dataset(options_.first_obs+options_.nobs:end, :);
        forecast_data.realized_xdata = ones(size(forecast_data.realized_val, 1), nx);
    else
        forecast_data.realized_val = [];
    end
end


function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
%function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
% ydum, xdum:   dummy observation data that implement the prior
% breaks:       vector of points in the dummy data after which new dummy obs's start
%                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where 
%                   actual data matrix has T rows, in preparing input for rfvar3
% nv,nx,lags: VAR dimensions
% mnprior.tight:Overall tightness of Minnesota prior
% mnprior.decay:Standard deviations of lags shrink as lag^(-decay)
% vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
% vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
%                   Should be an integer, and will be rounded if not.  vprior.sig is needed
%                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
%                   Set vprior.w=0 to achieve this.
% Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
%                   cannot be implemented entirely with dummy observations.  It is also usually
%                   taken to include the sum-of-coefficients and co-persistence components
%                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
%                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
%                   prior mean generates larger prior variances for own than for cross-effects even in 
%                   this formulation, but here there is no way to shrink toward a set of unconstrained 
%                   univariate AR's.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m

if ~isempty(mnprior)
    xdum = zeros(lags+1,nx,lags,nv);
    ydum = zeros(lags+1,nv,lags,nv);
    for il = 1:lags
        ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
    end
    ydum(1,:,1,:) = diag(vprior.sig);
    ydum = mnprior.tight*reshape(ydum,[lags+1,nv,lags*nv]);
    ydum = flipdim(ydum,1);
    xdum = mnprior.tight*reshape(xdum,[lags+1,nx,lags*nv]);
    xdum = flipdim(xdum,1);
    breaks = (lags+1)*[1:(nv*lags)]';
    lbreak = breaks(end);
else
    ydum = [];
    xdum = [];
    breaks = [];
    lbreak = 0;
end
if ~isempty(vprior) && vprior.w>0
    ydum2 = zeros(lags+1,nv,nv);
    xdum2 = zeros(lags+1,nx,nv);
    ydum2(end,:,:) = diag(vprior.sig);
    for i = 1:vprior.w
        ydum = cat(3,ydum,ydum2);
        xdum = cat(3,xdum,xdum2);
        breaks = [breaks;(lags+1)*[1:nv]'+lbreak];
        lbreak = breaks(end);
    end
end
dimy = size(ydum);
ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),nv);
xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
breaks = breaks(1:(end-1));
