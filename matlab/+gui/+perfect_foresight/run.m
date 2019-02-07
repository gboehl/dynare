function varargout = run(json)
% function varargout = run(json)
% Read JSON and run perfect foresight solver. Potentially return output as
% JSON
%
% INPUTS
%   json         [string]   JSON string representing options to run perfect
%                           foresight solver
%
% OUTPUTS
%   varargout{1} [string]   if desired, return output as JSON string
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_ oo_ ys0_ ex0_

%% Check Inputs
if nargin ~= 1 || ~ischar(json)
    error('function takes one string input argument')
end

if nargout > 1
    error('function provides up to one output argument')
end

%% Read JSON
jm = loadjson(json, 'SimplifyCell', 0);

%% INITVAL instructions
% initialize exogenous shocks to zero and compute initial steady state
options_.initval_file = 0;
oo_.exo_steady_state(:, 1) = 0;
if M_.exo_nbr > 0
    oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
    oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
if nargout == 1
    data2json = struct();
    data2json.steady_state1 = oo_.steady_state;
end

%% ENDVAL instructions
% initialize exogenous shocks to zero and compute final ss unless there is a permanent shock
M_.det_shocks = [];
ys0_= oo_.steady_state;
ex0_ = oo_.exo_steady_state;
permanent_shock_exists = isfield(jm, 'permanent_shocks') && ~isempty(jm.permanent_shocks);
if permanent_shock_exists
    for i = 1:length(jm.permanent_shocks)
        s = jm.permanent_shocks(i);
        oo_.exo_steady_state(s.index) = s.value;
        if s.start_period > 1
            % if the permanent shock does not start at the initial period
            % add a shocks block to mask the unnecessary periods
            M_.det_shocks = [ ...
                M_.det_shocks; ...
                struct(...
                'exo_det', 0, ...
                'exo_id', s.index, ...
                'multiplicative', 0, ...
                'periods', 1:s.start_period, ...
                'value', 0)];
        end
    end
else
    oo_.exo_steady_state(:, 1) = 0;
end
steady;
savedpermanentSS = oo_.steady_state;
if nargout == 1
    data2json.steady_state2 = oo_.steady_state;
end

%% SHOCKS instructions (for transitory shocks)
if jm.transitoryshockexist == 1
    for exotriter = 1:length(jm.shocksdescription)
        currenttrshock = jm.shocksdescription(exotriter);
        M_.det_shocks = [ M_.det_shocks;struct('exo_det',0,'exo_id',(currenttrshock{1}.shockindex+1),'multiplicative',0,'periods',currenttrshock{1}.shockstartperiod:currenttrshock{1}.shockendperiod,'value',currenttrshock{1}.shockvalue) ];
    end
    M_.exo_det_length = 0;
end

if jm.nonanticipatedshockexist == 1 || jm.delayexist == 1
    nonanticip = jm.nonanticipmatrix;
    rowindex = 1;
    firstsimul = 0;
    while nonanticip{rowindex}{1} > 0
        currentperiod=nonanticip{rowindex}{1};
        if currentperiod == 1
            % there are nonanticipated shocks to add at first period
            if nonanticip{rowindex}{4} == 0
                % this is a current nonanticipated shock
                M_.det_shocks = [ M_.det_shocks;struct('exo_det',0,'exo_id',(nonanticip{rowindex}{2}+1),'multiplicative',0,'periods',1:1,'value',nonanticip{rowindex}{7}) ];
            else
                % this is a delayed nonanticipated shock
                M_.det_shocks = [ M_.det_shocks;struct('exo_det',0,'exo_id',(nonanticip{rowindex}{2}+1),'multiplicative',0,'periods',(nonanticip{rowindex}{5}):(nonanticip{rowindex}{6}),'value',nonanticip{rowindex}{7}) ];
            end
            if nonanticip{rowindex+1}{1} ~= currentperiod
                % when we have tracked all first period shocks we can simulate
                options_.periods = jm.simperiods;
                yy = oo_.steady_state;
                perfect_foresight_setup;
                [rowexo, colexo] = size(oo_.exo_simul);
                perfect_foresight_solver;

                if nonanticip{rowindex+1}{1} > 0
                    % we collect all the path from ooendo period 1 to just before the next shock...
                    yy = [yy oo_.endo_simul(:,2:(2+(nonanticip{rowindex+1}{1}-currentperiod-1)))];
                else
                    % or if there are no more shocks we collect the whole path
                    yy = [yy oo_.endo_simul(:,2:end)];
                end

                ooexosaved = oo_.exo_simul;
                firstsimul = 1;
            end
        else
            % currentperiod is larger than one: we first perform perfect foresight simulation with initial period 1 conditions
            if firstsimul == 0
                % Initializing the first simulation
                options_.periods = jm.simperiods;
                yy = oo_.steady_state;
                perfect_foresight_setup;
                [rowexo, colexo] = size(oo_.exo_simul);
                perfect_foresight_solver;

                % In this because there is at least one shock we did not consider yet in the first period, we only save the path from the beginning up the period just before the current
                yy = [yy oo_.endo_simul(:,2:currentperiod)];
                ooexosaved = oo_.exo_simul;
                firstsimul = 1;
            end

            if nonanticip{rowindex}{3} == 1
                % permanent shock
                oo_.exo_steady_state(nonanticip{rowindex}{2}+1) = nonanticip{rowindex}{7};
                steady;
                savedpermanentSS = oo_.steady_state;
                if nargout == 1
                    data2json.steady_state2 = oo_.steady_state;
                end

                if nonanticip{rowindex}{4} == 0
                    % current permanent nonanticipated shock
                    ooexosaved(currentperiod+1:end, nonanticip{rowindex}{2}+1) = nonanticip{rowindex}{7};
                else
                    % delayed permanent nonanticipated shock
                    ooexosaved(nonanticip{rowindex}{5}+1:end, nonanticip{rowindex}{2}+1) = nonanticip{rowindex}{7};
                end

            else
                % not a permanent shock
                % add new shocks in the saved timepath with original time indexes
                if nonanticip{rowindex}{4} == 0
                    % this is a single current nonanticipated shock
                    ooexosaved(currentperiod+1, nonanticip{rowindex}{2}+1) = nonanticip{rowindex}{7};
                else
                    % this is a delayed nonanticipated shock
                    ooexosaved(nonanticip{rowindex}{5}+1:nonanticip{rowindex}{6}+1, nonanticip{rowindex}{2}+1) = nonanticip{rowindex}{7};
                end
            end

            % copy only the necessary window in oo_.exo_simul
            oo_.exo_simul = [zeros(1, colexo); ooexosaved(currentperiod+1:end, :)];

            % fill oo_.exo_simul until it has the correct size depending on of there are permanent shocks or not
            if permanent_shock_exists
                % if there is a permanent shock, fill with last value of ooexosaved
                oo_.exo_simul = [oo_.exo_simul; ones(rowexo-size(oo_.exo_simul, 1), 1)*ooexosaved(end, :)];
            else
                % otherwise fill with zeros
                oo_.exo_simul = [oo_.exo_simul; zeros(rowexo-size(oo_.exo_simul, 1), colexo)];
            end

            if nonanticip{rowindex+1}{1} ~= currentperiod
                % when we have tracked all the non-anticipated/delayed shocks for the current period, we can simulate
                if permanent_shock_exists
                    % if there are permanent shocks, fill oo_.endo with finalSS
                    oo_.endo_simul = savedpermanentSS*ones(1, options_.periods+2);
                else
                    % no permanent shocks, fill oo_.endo with initialSS
                    oo_.endo_simul = oo_.steady_state*ones(1, options_.periods+2);
                end

                % change oo_.endo_simul first value that gives the initial state of the economy
                oo_.endo_simul(:, 1) = yy(:,end);

                perfect_foresight_solver;
                if nonanticip{rowindex+1}{1} > 0
                    % collect all the path from ooendo period 1 to just before the next shock...
                    yy = [yy oo_.endo_simul(:, 2:2+nonanticip{rowindex+1}{1}-currentperiod-1)];
                else
                    % or if there are no more shocks we collect the whole path
                    yy = [yy oo_.endo_simul(:, 2:end)];
                end
            end
        end
        rowindex = rowindex+1;
    end
    % copy the endo path back
    oo_.endo_simul = yy;
else
    % if there are no unanticipated shocks we perform the simulation
    options_.periods = jm.simperiods;
    perfect_foresight_setup;
    perfect_foresight_solver;
end

if nargout == 1
    plotlgt = length(oo_.endo_simul);
    data2json.endosimul_length = plotlgt;
    data2json.endo_names = char(M_.endo_names);
    data2json.endo_nbr = M_.endo_nbr;
    for nendo = 1:M_.endo_nbr
        data2json.endo_simul.(strtrim(char(M_.endo_names(nendo, :)))) = oo_.endo_simul(nendo, :);
    end
    data2json.endo_simul.plotx = 0:plotlgt;
    varargout{1} = savejson('', data2json, '');
end
end
