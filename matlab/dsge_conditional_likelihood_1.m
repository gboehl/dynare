function [fval, info, exitflag, DLIK, Hess, SteadyState, trend_coeff, Model, DynareOptions, BayesInfo, DynareResults] = ...
    dsge_conditional_likelihood_1(xparam1, DynareDataset, DatasetInfo, DynareOptions, Model, EstimatedParameters, BayesInfo, BoundsInfo, DynareResults, derivatives_info)

% Copyright (C) 2017-2023 Dynare Team
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

% Initialization of the returned variables and others...
fval        = [];
SteadyState = [];
trend_coeff = [];
exitflag    = true;
info        = zeros(4,1);
DLIK        = [];
Hess        = [];

% Exit with error if analytical_derivation option is used.
if DynareOptions.analytic_derivation
    error('The analytic_derivation and conditional_likelihood are not compatible!')
end

% Ensure that xparam1 is a column vector.
% (Don't do the transformation if xparam1 is empty, otherwise it would become a
%  0×1 matrix, which create issues with older MATLABs when comparing with [] in
%  check_bounds_and_definiteness_estimation)
if ~isempty(xparam1)
    xparam1 = xparam1(:);
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------
Model = set_all_parameters(xparam1,EstimatedParameters,Model);

[fval, info, exitflag, Q, H] = check_bounds_and_definiteness_estimation(xparam1, Model, EstimatedParameters, BoundsInfo);
if info(1)
    return
end

iQ_upper_chol = chol(inv(Q));

% Return an error if the interface for measurement errors is used.
if ~isequal(H, zeros(size(H))) || EstimatedParameters.ncn || EstimatedParameters.ncx
    error('Option conditional_likelihood does not support declaration of measurement errors. You can specify the measurement errors in the model block directly by adding measurement equations.')
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic steadystate and extract the matrices of the state equation (T and R).
[T, R, SteadyState, info,DynareResults.dr, Model.params] = ...
    dynare_resolve(Model, DynareOptions, DynareResults.dr, DynareResults.steady_state, DynareResults.exo_steady_state, DynareResults.exo_det_steady_state, 'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86 || ...
                info(1) == 401 || info(1) == 402 || info(1) == 403 || ... %cycle reduction
                info(1) == 411 || info(1) == 412 || info(1) == 413 % logarithmic reduction 
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exitflag = false;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exitflag = false;
        return
    end
end

% check endogenous prior restrictions
info = endogenous_prior_restrictions(T, R, Model, DynareOptions, DynareResults);
if info(1)
    fval = Inf;
    info(4)=info(2);
    exitflag = false;
    return
end

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the constant vector of the measurement equation.
if ~DynareOptions.noconstant
    if DynareOptions.loglinear
        constant = log(SteadyState(BayesInfo.mfys));
    else
        constant = SteadyState(BayesInfo.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    [trend_addition, trend_coeff] = compute_trend_coefficients(Model, DynareOptions, DynareDataset.vobs, DynareDataset.nobs);
    Y = bsxfun(@minus, transpose(DynareDataset.data), constant)-trend_addition;
else
    trend_coeff = zeros(DynareDataset.vobs, 1);
    if ~DynareOptions.noconstant
        Y = bsxfun(@minus, transpose(DynareDataset.data), constant);
    else
        Y = transpose(DynareDataset.data);
    end
end

% Return an error if some observations are missing.
if DatasetInfo.missing.state
    error('Option conditional_likelihood is not compatible with missing observations.')
end

% Get the selection matrix (vector of row indices for T and R)
Z = BayesInfo.mf;

% Get the number of observed variables.
pp = DynareDataset.vobs;

% Get the number of variables in the state equations (state variables plus observed variables).
mm = size(T, 1);

% Get the number of innovations.
rr = length(Q);

% Return an error if the number of shocks is not equal to the number of observations.
if ~isequal(pp, rr)
    error('With conditional_likelihood the number of innovations must be equal to the number of observed varilables!')
end

% Set state vector (deviation to steady state)
S = zeros(mm, 1);

%------------------------------------------------------------------------------
% 3. Evaluate the conditional likelihood
%------------------------------------------------------------------------------

[L, U] = lu(R(Z,:)); % note that det(L)={-1,1} depending on the number of permutations so we can forget it when we take the absolute value of the determinant of R(Z,:) below (in the constant).

const = -.5*rr*log(2*pi) - log(abs(prod(diag(U)))) + sum(log(diag(iQ_upper_chol)));

llik = zeros(size(Y, 2), 1);

Ytild = U\(L\Y);
Ttild = U\(L\T(Z,:));

for t = 1:DynareOptions.presample
    epsilon = Ytild(:,t) - Ttild*S;
    S = T*S + R*epsilon;
end

for t=(DynareOptions.presample+1):size(Y, 2)
    epsilon = Ytild(:,t) - Ttild*S;
    upsilon = iQ_upper_chol*epsilon;
    S = T*S + R*epsilon;
    llik(t) = const - .5*dot(upsilon, upsilon);
end

% Computes minus log-likelihood.
likelihood = -sum(llik);


% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------

lnprior = priordens(xparam1, BayesInfo.pshape, BayesInfo.p6, BayesInfo.p7, BayesInfo.p3, BayesInfo.p4);

if DynareOptions.endogenous_prior==1
    [lnpriormom]  = endogenous_prior(Y, Pstar, BayesInfo, H);
    fval = (likelihood-lnprior-lnpriormom);
else
    fval = (likelihood-lnprior);
end

if DynareOptions.prior_restrictions.status
    tmp = feval(DynareOptions.prior_restrictions.routine, Model, DynareResults, DynareOptions, DynareDataset, DatasetInfo);
    fval = fval - tmp;
end

if isnan(fval)
    fval = Inf;
    info(1) = 47;
    info(4) = 0.1;
    exitflag = false;
    return
end

if imag(fval)~=0
    fval = Inf;
    info(1) = 48;
    info(4) = 0.1;
    exitflag = false;
    return
end