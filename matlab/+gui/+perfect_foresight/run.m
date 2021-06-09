function run(json)
% function varargout = run(json)
% Read JSON and run perfect foresight solver. Potentially return output as
% JSON
%
% INPUTS
%   json         [string]   JSON string representing options to run perfect
%                           foresight solver
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ options_ oo_ ys0_ ex0_

%% Check Inputs
if nargin ~= 1 || ~ischar(json)
    error('function takes one string input argument')
end

if nargout > 1
    error('function provides up to one output argument')
end

%% Read JSON
jm = loadjson_(json, 'SimplifyCell', 1);

%% INITVAL instructions
% initialize exogenous shocks to zero and compute initial steady state
options_.initval_file = 0;
oo_.steady_state(:, 1) = 0;
for i = 1:length(jm.initval_endo)
    oo_.steady_state(jm.initval_endo(i).id) = jm.initval_endo(i).value;
end
oo_.exo_steady_state(:, 1) = 0;
for i = 1:length(jm.initval_exo)
    oo_.exo_steady_state(jm.initval_exo(i).id) = jm.initval_exo(i).value;
end
if M_.exo_nbr > 0
    oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
    oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end

%% ENDVAL instructions
% initialize exogenous shocks to zero and compute final ss unless there is a permanent shock
ys0_ = [];
ex0_ = [];
M_.det_shocks = [];
if ~isempty(jm.anticipated_permanent_shocks) || ~isempty(jm.endval_endo)
    ys0_= oo_.steady_state;
    ex0_ = oo_.exo_steady_state;
    for i = 1:length(jm.endval_endo)
        oo_.steady_state(jm.endval_endo(i).id) = jm.endval_endo(i).value;
    end
    for i = 1:length(jm.anticipated_permanent_shocks)
        s = jm.anticipated_permanent_shocks(i);
        oo_.exo_steady_state(s.exo_id) = s.value;
        if s.start_date > 1
            % if the permanent shock does not start at the initial period
            % add a shocks block to mask the unnecessary periods
            M_.det_shocks = [ ...
                M_.det_shocks; ...
                struct(...
                'exo_det', 0, ...
                'exo_id', s.exo_id, ...
                'multiplicative', 0, ...
                'periods', 1:s.start_date, ...
                'value', 0)];
        end
    end
end

%% SHOCKS instructions (for anticipated transitory shocks)
if ~isempty(jm.anticipated_transitory_shocks)
    for i = 1:length(jm.anticipated_transitory_shocks)
        s = jm.anticipated_transitory_shocks(i);
        M_.det_shocks = [ ...
            M_.det_shocks; ...
            struct('exo_det', 0, ...
            'exo_id', s.exo_id, ...
            'multiplicative', 0, ...
            'periods', s.start_date:s.end_date, ...
            'value', s.value)];
    end
    M_.exo_det_length = 0;
end

%% Make unanticipated shock map
unanticipated_p_shocks = containers.Map('KeyType', 'int32', 'ValueType', 'any');
for i = 1:length(jm.unanticipated_permanent_shocks)
    s = jm.unanticipated_permanent_shocks(i);
    if isempty(s.anticipated_date)
        unanticipated_p_shocks(s.start_date) = s;
    else
        if s.anticipated_date > s.start_date
            error('The expected date cannot be greater than the shock start date')
        end
        unanticipated_p_shocks(s.anticipated_date) = s;
    end
end

unanticipated_t_shocks = containers.Map('KeyType', 'int32', 'ValueType', 'any');
for i = 1:length(jm.unanticipated_transitory_shocks)
    s = jm.unanticipated_transitory_shocks(i);
    if isempty(s.anticipated_date)
        for j = s.start_date:s.end_date
            ts = s;
            ts.start_date = j;
            ts.end_date = j;
            unanticipated_t_shocks(j) = ts;
        end
    else
        if s.anticipated_date > s.start_date
            error('The expected date cannot be greater than the shock start date')
        end
        unanticipated_t_shocks(s.anticipated_date) = s;
    end
end

mapkeys = unique(cell2mat([keys(unanticipated_p_shocks) keys(unanticipated_t_shocks)]));

%% Simulation
options_.periods = jm.periods;
perfect_foresight_setup;

% no surprise shocks present
if isempty(mapkeys)
    perfect_foresight_solver;
    return
end

% surprise shocks present
% in case there are unanticipated shocks...
if isempty(ys0_)
    yy = oo_.steady_state;
else
    yy = ys0_;
end

if mapkeys(1) ~= 1
    % if first unanticipated shock is not in period 1
    % simulate until first unanticipated shock and save
    perfect_foresight_solver;
    yy = [yy oo_.endo_simul(:, 2:mapkeys(1)+1)];
end

last_period = 1;
length(oo_.exo_simul)
oo_exo_simul_rows = options_.periods + 2;
for i = 1:length(mapkeys)
    this_period = mapkeys(i);
    if i ~= length(mapkeys)
        next_period = mapkeys(i+1);
    else
        next_period = -1;
    end
    if mapkeys(i) ~= 1
        % shift shock path
        nperiods = this_period - last_period;
        oo_.exo_simul = [oo_.exo_simul(nperiods+1:end, :); repmat(oo_.exo_steady_state, nperiods, 1)];
    end
    if isKey(unanticipated_p_shocks, mapkeys(i))
        s = unanticipated_p_shocks(mapkeys(i));
        if isempty(s.anticipated_date) || s.start_date == s.anticipated_date
            oo_.exo_steady_state(s.exo_id) = s.value;
            oo_.exo_simul(2:end, :) = repmat(oo_.exo_steady_state, oo_exo_simul_rows-1, 1);
        else
            date_offset = s.start_date - s.anticipated_date;
            oo_.exo_steady_state(s.exo_id) = s.value;
            oo_.exo_simul(date_offset+1:end, :) = repmat(oo_.exo_steady_state, oo_exo_simul_rows-date_offset-1, 1);
        end
    end
    if isKey(unanticipated_t_shocks, mapkeys(i))
        s = unanticipated_t_shocks(mapkeys(i));
        if isempty(s.anticipated_date) || s.start_date == s.anticipated_date
            oo_.exo_simul(2, s.exo_id) = s.value;
        else
            date_offset = s.start_date - s.anticipated_date;
            oo_.exo_simul(date_offset+1:s.end_date-s.start_date+1+date_offset, s.exo_id) = s.value;
        end
    end
    last_period = this_period;
    assert(rows(oo_.exo_simul) == oo_exo_simul_rows, 'error encountered setting oo_.exo_simul');
    oo_.endo_simul(:, 1) = yy(:, end);
    perfect_foresight_solver;
    if next_period > 0
        yy = [yy oo_.endo_simul(:, 2:next_period-this_period+1)];
    else
        assert(i == length(mapkeys), 'should not arrive here');
        yy = [yy oo_.endo_simul(:, 2:end)];
    end
end
oo_.endo_simul = yy;
end
