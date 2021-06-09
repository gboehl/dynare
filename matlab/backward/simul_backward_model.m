function simulation = simul_backward_model(initialconditions, samplesize, innovations)

% Simulates a stochastic backward looking model (with arbitrary precision).
%
% INPUTS
% - initialconditions   [double]      n*1 vector, initial conditions for the endogenous variables.
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
% [3] If the first input argument is empty, the endogenous variables are initialized with 0, or if available with the informations
%     provided thrtough the histval block.

% Copyright (C) 2012-2019 Dynare Team
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

global options_ M_ oo_

if M_.maximum_lead
    error('Model defined in %s.mod is not backward or static.', M_.fname)
end

if ~M_.maximum_lag
    dprintf('Model defined in %s.mod is static. Use simul_static_model instead.', M_.fname)
    simul_static_model(samplesize, innovations);
    return
end

if ismember(options_.solve_algo, [12,14]) && ~M_.possible_to_use_solve_algo_12_14
    error(M_.message_solve_algo_12_14)
end

if nargin<3
    Innovations = [];
else
    if isdseries(innovations)
        if isdseries(initialconditions)
            if isequal(innovations.dates(1)-1, initialconditions.dates(end))
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
                        dprintf('Exogenous variable %s is not available in third argument, default value is zero.', exonames{i})
                    end
                end
            else
                error('Time spans in first and third arguments should be contiguous!')
            end
        else
            if isempty(initialconditions)
                if innovations.nobs<samplesize
                    error('Time span in third argument is too short (should not be less than %s, the value of the second argument)', num2str(samplesize))
                end
                Innovations = zeros(samplesize, M_.exo_nbr);
                exonames = M_.exo_names;
                for i=1:M_.exo_nbr
                    if ismember(exonames{i}, innovations.name)
                        Innovations(:,i) = innovations{exonames{i}}.data(1:samplesize);
                    else
                        dprintf('Exogenous variable %s is not available in third argument, default value is zero.', exonames{i})
                    end
                end
            else
                error('First input must be an empty array!')
            end
        end
    else
        error('Third argument must be a dseries object!')
    end
end

if options_.linear
    simulation = simul_backward_linear_model(initialconditions, samplesize, options_, M_, oo_, Innovations);
else
    simulation = simul_backward_nonlinear_model(initialconditions, samplesize, options_, M_, oo_, Innovations);
end