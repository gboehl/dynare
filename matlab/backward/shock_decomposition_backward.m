function decomposition = shock_decomposition_backward(simulations, initialconditions, shocklist, endograph, use_shock_groups)

% Computes and possibly plots the shock decomposition of a backward (possibly nonlinear)
% model simulation.
%
% Inputs:
% - simulations          dseries object as returned by simul_backward_model
% - initialconditions    dseries object as passed to simul_backward_model
% - shocklist            a cell array listing the (names of the) shocks whose contribution should be
%                        computed. The order matters: the contribution of a shock appearing at index
%                        i is computed as the difference between the simulation where all shocks ≥i+1
%                        are zero and the simulation where all shocks ≥i are zero. It is also
%                        possible to put in this list shock groups, as defined in a shock_groups
%                        block
% - endograph            an optional cell array listing the (names of the) endogenous variables for
%                        which a shock decomposition graph should be created
% - use_shock_groups     an optional string giving the name of the shock_groups block to be used to
%                        resolve shock groups listed in shocklist. If not given, 'default' is used.
%
% Output:
% - decomposition        a 3D array of size endo_nbr×nshocks×nperiods where nshocks=length(shocklist)
%                        and nperiods is the number of simulation periods (i.e. excluding the initial
%                        conditions)

% Copyright © 2020-2021 Dynare Team
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

global M_ options_

narginchk(3, 5);
if nargin < 4
    endograph = {};
end
if nargin < 5
    use_shock_groups = 'default';
end

% Extract matrix of innovations from the whole simulation paths
simdates = simulations.dates(initialconditions.nobs+1:end);
innovations = simulations{M_.exo_names{:}}(simdates);

% Number of simulation periods
nperiods = simulations.nobs - initialconditions.nobs;

decomposition = NaN(M_.endo_nbr, length(shocklist), nperiods);

% Add auxiliary variables to simulation paths
if exist(sprintf('+%s/dynamic_set_auxiliary_series.m', M_.fname), 'file')
    simulations = feval(sprintf('%s.dynamic_set_auxiliary_series', M_.fname), simulations, M_.params);
end

for i = length(shocklist):-1:1
    % Zero the innovation(s) corresponding to this shock or shock group
    if ismember(shocklist{i}, M_.exo_names)
        innovations{shocklist{i}}(simdates) = zeros(nperiods, 1);
    else % This is a shock group
        if ~ismember(use_shock_groups, fieldnames(M_.shock_groups))
            error(['Unknown shock_groups block: ' use_shock_groups])
        end
        groups = fieldnames(M_.shock_groups.(use_shock_groups));
        shocks_in_group = [];
        for j = 1:length(groups)
            if strcmp(shocklist{i}, M_.shock_groups.(use_shock_groups).(groups{j}).label)
                shocks_in_group = M_.shock_groups.(use_shock_groups).(groups{j}).shocks;
                break
            end
        end
        if isempty(shocks_in_group)
            error(['Unknown shock group: ' shocklist{i}])
        end
        for j = 1:length(shocks_in_group)
            innovations{shocks_in_group{j}}(simdates) = zeros(nperiods, 1);
        end
    end

    % Compute simulation with the current shock or shock group removed
    simulations_new = simul_backward_model(initialconditions, nperiods, innovations);
    if exist(sprintf('+%s/dynamic_set_auxiliary_series.m', M_.fname), 'file')
        simulations_new = feval(sprintf('%s.dynamic_set_auxiliary_series', M_.fname), simulations_new, M_.params);
    end

    % Compute the contribution of the current shock or shock group
    if ~isoctave && matlab_ver_less_than('9.7')
        % Workaround for MATLAB bug described in dseries#45
        % The solution is to avoid using the "end" keyword
        myend = nobs(simulations);
        contribution = simulations{M_.endo_names{:}}.data(initialconditions.nobs+1:myend, :) ...
                       - simulations_new{M_.endo_names{:}}.data(initialconditions.nobs+1:myend, :);
    else
        contribution = simulations{M_.endo_names{:}}.data(initialconditions.nobs+1:end, :) ...
                       - simulations_new{M_.endo_names{:}}.data(initialconditions.nobs+1:end, :);
    end
    decomposition(:, i, :) = contribution';

    simulations = simulations_new;
end

% Plot the decomposition

for i = 1:length(endograph)
    h = dyn_figure(options_.plot_shock_decomp.nodisplay, 'Name', [ 'Shock decomposition for ' endograph{i}]);
    endoidx = find(strcmp(endograph{i}, M_.endo_names));
    bar(double(simdates), squeeze(decomposition(endoidx, :, :))', 'stacked')
    legend(shocklist)
end

end
