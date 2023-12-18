function Simulations = extended_path_mc(initialconditions, samplesize, replic, exogenousvariables, options_, M_, oo_)
% Simulations = extended_path_mc(initialconditions, samplesize, replic, exogenousvariables, options_, M_, oo_)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
%  o initialconditions      [double]    m*1 array, where m is the number of endogenous variables in the model.
%  o samplesize             [integer]   scalar, size of the sample to be simulated.
%  o exogenousvariables     [double]    T*n array, values for the structural innovations.
%  o options_               [struct]    Dynare's options structure
%  o M_                     [struct]    Dynare's model structure
%  o oo_                    [struct]    Dynare's results structure
%
% OUTPUTS
%  o ts                     [dseries]   m*samplesize array, the simulations.
%  o results                [cell]
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright © 2016-2023 Dynare Team
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

[initialconditions, innovations, pfm, ep, ~, options_, oo_] = ...
    extended_path_initialization(initialconditions, samplesize, exogenousvariables, options_, M_, oo_);

% Check the dimension of the first input argument
if isequal(size(initialconditions, 2), 1)
    initialconditions = repmat(initialconditions, 1, replic);
else
    if ~isequal(size(initialconditions, 2), replic)
        error('Wrong size. Number of columns in first argument should match the value of the third argument!')
    end
end

% Check the dimension of the fourth input argument
if isempty(exogenousvariables)
    exogenousvariables = repmat(exogenousvariables, [1, 1, replic]);
else
    if ~isequal(size(exogenousvariables, 3), replic)
        error('Wrong size. !')
    end
end
if ~isequal(size(exogenousvariables, 3), replic)
    error('Wrong dimensions. Fourth argument must be a 3D array with as many pages as the value of the third argument!')
end

data = NaN(size(initialconditions, 1), samplesize+1, replic);
vexo = NaN(innovations.effective_number_of_shocks, samplesize+1, replic);
info = NaN(replic, 1);

if ep.parallel
    % Use the Parallel toolbox.
    parfor i=1:replic
        innovations_ = innovations;
        oo__ = oo_;
        [shocks, spfm_exo_simul, innovations_, oo__] = extended_path_shocks(innovations_, ep, exogenousvariables(:,:,i), samplesize, M_, options_, oo__);
        endogenous_variables_paths = NaN(M_.endo_nbr,samplesize+1);
        endogenous_variables_paths(:,1) = initialconditions(:,1);
        exogenous_variables_paths = NaN(innovations_.effective_number_of_shocks,samplesize+1);
        exogenous_variables_paths(:,1) = 0;
        info_convergence = true;
        t = 1;
        while t<=samplesize
            t = t+1;
            spfm_exo_simul(2,:) = shocks(t-1,:);
            exogenous_variables_paths(:,t) = shocks(t-1,:);
            [endogenous_variables_paths(:,t), info_convergence] = extended_path_core(ep.periods, M_.endo_nbr, M_.exo_nbr, innovations_.positive_var_indx, ...
                                                              spfm_exo_simul, ep.init, endogenous_variables_paths(:,t-1), ...
                                                              oo__.steady_state, ...
                                                              ep.verbosity, ep.stochastic.order, ...
                                                              M_, pfm,ep.stochastic.algo, ep.solve_algo, ep.stack_solve_algo, ...
                                                              options_.lmmcp, options_, oo__);
            if ~info_convergence
                msg = sprintf('No convergence of the (stochastic) perfect foresight solver (in period %s, iteration %s)!', int2str(t), int2str(i));
                warning(msg)
                break
            end
        end % Loop over t
        info(i) = info_convergence;
        vexo(:,:,i) = exogenous_variables_paths;
        data(:,:,i) = endogenous_variables_paths;
    end
else
    % Sequential approach.
    for i=1:replic
        [shocks, spfm_exo_simul, innovations, oo_] = extended_path_shocks(innovations, ep, exogenousvariables(:,:,i), samplesize, M_, options_, oo_);
        endogenous_variables_paths = NaN(M_.endo_nbr,samplesize+1);
        endogenous_variables_paths(:,1) = initialconditions(:,1);
        exogenous_variables_paths = NaN(innovations.effective_number_of_shocks,samplesize+1);
        exogenous_variables_paths(:,1) = 0;
        t = 1;
        while t<=samplesize
            t = t+1;
            spfm_exo_simul(2,:) = shocks(t-1,:);
            exogenous_variables_paths(:,t) = shocks(t-1,:);
            [endogenous_variables_paths(:,t), info_convergence] = extended_path_core(ep.periods, M_.endo_nbr, M_.exo_nbr, innovations.positive_var_indx, ...
                                                              spfm_exo_simul, ep.init, endogenous_variables_paths(:,t-1), ...
                                                              oo_.steady_state, ...
                                                              ep.verbosity, ep.stochastic.order, ...
                                                              M_, pfm,ep.stochastic.algo, ep.solve_algo, ep.stack_solve_algo, ...
                                                              options_.lmmcp, options_, oo_);
            if ~info_convergence
                msg = sprintf('No convergence of the (stochastic) perfect foresight solver (in period %s, iteration %s)!', int2str(t), int2str(i));
                warning(msg)
                break
            end
        end % Loop over t
        info(i) = info_convergence;
        vexo(:,:,i) = exogenous_variables_paths;
        data(:,:,i) = endogenous_variables_paths;
    end % Loop over i
end

Simulations.innovations = vexo;
Simulations.data = data;
Simulations.info = info;
