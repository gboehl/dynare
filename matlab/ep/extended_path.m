function [ts,oo_] = extended_path(initialconditions, samplesize, exogenousvariables, options_, M_, oo_)
% [ts,oo_] = extended_path(initialconditions, samplesize, exogenousvariables, options_, M_, oo_)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
%  o initialconditions      [double]    m*1 array, where m is the number of endogenous variables in the model.
%  o samplesize             [integer]   scalar, size of the sample to be simulated.
%  o exogenousvariables     [double]    T*n array, values for the structural innovations.
%  o options_               [struct]    options_
%  o M_                     [struct]    Dynare's model structure
%  o oo_                    [struct]    Dynare's results structure
%
% OUTPUTS
%  o ts                     [dseries]   m*samplesize array, the simulations.
%  o results                [struct]    results structure
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright © 2009-2023 Dynare Team
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

[initialconditions, innovations, pfm, ep, verbosity, options_, oo_] = ...
    extended_path_initialization(initialconditions, samplesize, exogenousvariables, options_, M_, oo_);

[shocks, spfm_exo_simul, innovations, oo_] = extended_path_shocks(innovations, ep, exogenousvariables, samplesize,M_,options_,oo_);

% Initialize the matrix for the paths of the endogenous variables.
endogenous_variables_paths = NaN(M_.endo_nbr,samplesize+1);
endogenous_variables_paths(:,1) = initialconditions;

% Set waitbar (graphic or text  mode)
hh_fig = dyn_waitbar(0,'Please wait. Extended Path simulations...');
set(hh_fig,'Name','EP simulations.');

% Initialize while-loop index.
t = 1;

% Main loop.
while (t <= samplesize)
    if ~mod(t,10)
        dyn_waitbar(t/samplesize,hh_fig,'Please wait. Extended Path simulations...');
    end
    % Set period index.
    t = t+1;
    spfm_exo_simul(2,:) = shocks(t-1,:);
    if t>2
        % Set initial guess for the solver (using the solution of the
        % previous period problem).
        initialguess = [endogenousvariablespaths(:, 2:end), oo_.steady_state];
    else
        initialguess = [];
    end
    [endogenous_variables_paths(:,t), info_convergence, endogenousvariablespaths] = extended_path_core(ep.periods, M_.endo_nbr, M_.exo_nbr, innovations.positive_var_indx, ...
                                                      spfm_exo_simul, ep.init, endogenous_variables_paths(:,t-1), ...
                                                      oo_.steady_state, ...
                                                      verbosity, ep.stochastic.order, ...
                                                      M_, pfm, ep.stochastic.algo, ep.solve_algo, ep.stack_solve_algo, ...
                                                      options_.lmmcp, ...
                                                      options_, ...
                                                      oo_, initialguess);
    if ~info_convergence
        msg = sprintf('No convergence of the (stochastic) perfect foresight solver (in period %s)!', int2str(t));
        warning(msg)
        break
    end
end % (while) loop over t

% Close waitbar.
dyn_waitbar_close(hh_fig);

% Set the initial period.
if isdates(options_.initial_period)
    if ischar(options_.initial_period)
        initial_period = dates(options_.initial_period);
    else
        initial_period = options_.initial_period;
    end
elseif isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    error('Type of option initial_period is wrong.')
end

% Return the simulated time series.
if any(isnan(endogenous_variables_paths(:)))
    sl = find(~isnan(endogenous_variables_paths));
    nn = size(endogenous_variables_paths, 1);
    endogenous_variables_paths = reshape(endogenous_variables_paths(sl), nn, length(sl)/nn);
end
ts = dseries(transpose(endogenous_variables_paths), initial_period, M_.endo_names);

oo_.endo_simul = transpose(ts.data);