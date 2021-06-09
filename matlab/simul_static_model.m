function simulation = simul_static_model(samplesize, innovations)

% Simulates a stochastic static model (with arbitrary precision).
%
% INPUTS
% - samplesize          [integer]     scalar, number of periods for the simulation.
% - innovations         [dseries]     innovations to be used for the simulation.
%
% OUTPUTS
% - simulation          [dseries]     Simulated endogenous and exogenous variables.
%
% REMARKS
% [1] The innovations used for the simulation are saved in DynareOutput.exo_simul, and the resulting paths for the endogenous
%     variables are saved in DynareOutput.endo_simul.
% [2] The last input argument is not mandatory. If absent we use random draws and rescale them with the informations provided
%     through the shocks block.

% Copyright (C) 2019 Dynare Team
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

global M_ options_ oo_

if M_.maximum_lag
    error('%s.mod has lagged variables, but it should be a static model.', M_.fname)
end

if M_.maximum_lead
    error('%s.mod has leaded variables, but it should be a static model.', M_.fname)
end

% Set innovations.
if nargin<2 || isempty(innovations)
    % Set the covariance matrix of the structural innovations.
    variances = diag(M_.Sigma_e);
    number_of_shocks = length(M_.Sigma_e);
    positive_var_indx = find(variances>0);
    effective_number_of_shocks = length(positive_var_indx);
    covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx);
    covariance_matrix_upper_cholesky = chol(covariance_matrix);
    % Set seed to its default state.
    if options_.bnlms.set_dynare_seed_to_default
        set_dynare_seed('default');
    end
    % Simulate structural innovations.
    switch options_.bnlms.innovation_distribution
      case 'gaussian'
        oo_.bnlms.shocks = randn(samplesize, effective_number_of_shocks)*covariance_matrix_upper_cholesky;
      otherwise
        error('%s distribution for the structural innovations is not (yet) implemented!', options_.bnlms.innovation_distribution)
    end
    % Put the simulated innovations in DynareOutput.exo_simul.
    oo_.exo_simul = zeros(samplesize, number_of_shocks);
    oo_.exo_simul(:,positive_var_indx) = oo_.bnlms.shocks;
    innovations = [];
else
    if innovations.nobs<samplesize
        error('Time span in third argument is too short (should not be less than %s, the value of the second argument)', num2str(samplesize))
    end
    % Set array holding innovations values.
    Innovations = zeros(samplesize, M_.exo_nbr);
    exonames = M_.exo_names;
    for i=1:M_.exo_nbr
        if ismember(exonames{i}, innovations.name)
            Innovations(:,i) = innovations{exonames{i}}.data(1:samplesize);
        else
            dprintf('Exogenous variable %s is not available in third argument, default value is zero.', exonames{i});
        end
    end
    oo_.exo_simul = Innovations;
end

staticmodel = str2fun(sprintf('%s.static', M_.fname));

% Simulations (call a Newton-like algorithm for each period).
for t=1:samplesize
    y = zeros(M_.endo_nbr, 1);
    [oo_.endo_simul(:,t), info] = dynare_solve(staticmodel, y, options_, oo_.exo_simul(t,:), M_.params);
    if info
        error('Newton failed!')
    end
end

ysim = oo_.endo_simul(1:M_.orig_endo_nbr,:);
xsim = oo_.exo_simul;

initperiod = dates('1Y');
if isdseries(innovations)
    initperiod = innovations.dates(1);
end

simulation = [dseries(ysim', initperiod, M_.endo_names(1:M_.orig_endo_nbr)), dseries(xsim, initperiod, M_.exo_names)];
