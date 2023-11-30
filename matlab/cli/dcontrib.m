function dcontrib(varargin)

% Computes dynamic contributions to a subset of endogenous variables in a semi structural model.
%
% EXAMPLE
%
% >> dcontrib --model sandbox.mod --tags zpac eq:x1 --database ds --output results --range 2023Q1:2073Q1
%
% zpac and eq:x1 are the equation tags of the equations determining the endogenous variables for which we want to compute
% the contributions of the other (exogenous) variables, sandbox.mod is the name of the file from which we exctract these
% equations, ds is a dseries object containing the data, 2023Q1:2073Q1 is the time range over which we compute the
% contributions, and results the name of the structure containing the contributions (as dseries objects) for each endogenous
% variable.
%
% INPUTS
% --model                name of a mod file (with extension)
% --tags                 list of equations (equation tags assocated to the endogenous variables for which we want to compute the contributions)
% --database             dseries object
% --baseline             dseries object (path for the exogenous variables)
% --range                followed by a dates range
% --method               followed by cumulate (default) or diff.
% --log                  returns the variables in logs
%
% REMARKS
% [1] --baseline and --range are not compatible.
% [2] --variables is followed by a space separated list of names, it is assumed that each variable is associated with an equation tag.
% [3] In the context of an error correction or PAC equation, if one is willing to decompose the endogenous variable and the target then the equation tag for the target must also be provided.

% Copyright © 2023 Dynare Team
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

    global M_

    if nargin==1 && strcmpi(varargin{1}, '--help')
        skipline()
        disp('--model     followed by the name of a mod file (with extension) [mandatory]')
        disp('--tags      followed by a list of equation tags [mandatory]')
        disp('--database  followed by dseries object [mandatory]')
        disp('--baseline  followed by dseries object (path for the exogenous variables)')
        disp('--range     followed by a dates range')
        disp('--method    followed by keywords cumulate or diff')
        disp('--output    followed by a name for the structure holding the results [mandatory]')
        disp('--log       to return the contributions in logs')
        skipline()
        return
    end

    model = getmodel(varargin);

    % First call to dynare to obtain the json verison of the model.
    dynare(model(1:end-4), 'nopreprocessoroutput', 'notime', 'json=compute')
    delete(sprintf('%s.log', model(1:end-4)))

    eqtags = geteqtags(varargin);
    variables = cell(length(eqtags), 1);
    for i=1:length(eqtags)
        variables(i) = get_variables_and_parameters_in_expression(get_lhs_and_rhs(eqtags{i}, M_, true));
    end

    % Cherry pick equations required for the decomposition.
    cherrypickdir = sprintf('cherry-pick-%s', randomstring(10));
    cherrypick(model(1:end-4), cherrypickdir, eqtags, false);
    rmdir(model(1:end-4), 's')
    rmdir(sprintf('+%s', model(1:end-4)), 's')
    modfilename = sprintf('dcontrib_%s.mod', randomstring(10));
    aggregate(modfilename, {}, '', cherrypickdir);
    rmdir(cherrypickdir, 's')

    % Second call to dynare (on the exctracted equations)
    dynare(modfilename(1:end-4), 'nopreprocessoroutput', 'notime', 'json=compute')

    % Get dataset
    dname = getdatasetname(varargin);
    ds = evalin('caller', dname);
    if ~isdseries(ds)
        error('dcontrib:getdataset: --dataset must be followed by a dseries object.')
    end

    % Create a dseries object for the paths of the exogenous variables
    xvariables = ds{M_.exo_names{:}};

    % Get initial and terminal periods (if defined)
    [firstperiod, lastperiod] = getperiods(varargin);
    if firstperiod<=ds.dates(1)+M_.orig_maximum_lag
        error('dcontrib:: Try increase firstperiod (>%s).', char(ds.dates(1)+M_.orig_maximum_lag))
    end
    if lastperiod>ds.dates(end)
        error('dcontrib:: Try reduce lastperiod (<=%s).', char(ds.dates(end)))
    end

    if islog(varargin)
        transform = @(x) log(x);
    else
        transform = @(x) x;
    end

    % Load baseline (if it makes sense)
    if isempty(firstperiod)
        baselinename = getbaselinename(varargin);
        baseline = evalin('caller', baselinename);
        if ~isdseries(baseline)
            error('dcontrib:getdataset: --baseline must be followed by a dseries object.')
        end
        firstperiod = baseline.dates(1);
        lastperiod = baseline.dates(end);
        baseline = baseline{M_.exo_names{:}};
    else
        % Set default baseline (exogenous variable levels in firstperiod)
        baseline = xvariables(firstperiod);
        baseline = repmat(baseline.data, lastperiod-firstperiod+1, 1);
        baseline = dseries(baseline, firstperiod, M_.exo_names);
    end

    % get method for computing contributions
    method = getmethod(varargin);

    % Restrict the observations for the exogenous variables to the pertinent tim range
    xvariables = xvariables(firstperiod:lastperiod);

    % Set initial conditions for the simulation.
    initialconditions = ds(ds.dates(1):firstperiod-1);

    % Simulation on the baseline (track the effects of the initial state if the model is autoregressive)
    S.baseline = simul_backward_model(initialconditions, lastperiod-firstperiod+1, baseline);

    % contributions is a dseries object holding the marginal contribution of the baseline and
    % each exogenous variable to endogenous variable z

    switch method
      case 'cumulate'
        % Add exogenous variables one by one and simulate the model (-> cumulated contributions)
        for i=1:xvariables.vobs
            name = xvariables.name{i};
            baseline{name} = xvariables{name};
            S.(name) = simul_backward_model(initialconditions, lastperiod-firstperiod+1, baseline);
        end
        % Compute marginal contributions
        for j=1:length(variables)
            cumulatedcontribs = S.baseline{variables{j}}(firstperiod:lastperiod).data;
            contributions.(variables{j}) = dseries(transform(cumulatedcontribs), firstperiod, 'baseline');
            for i=1:xvariables.vobs
                name = xvariables.name{i};
                ts = S.(name);
                data = ts{variables{j}}(firstperiod:lastperiod).data;
                contributions.(variables{j}) = [contributions.(variables{j}), dseries(transform(data)-transform(cumulatedcontribs), firstperiod, name)];
                cumulatedcontribs = data;
            end
            contributions.(variables{j}) = contributions.(variables{j})(firstperiod:lastperiod);
        end
      case 'diff'
        for i=1:xvariables.vobs
            name = xvariables.name{i};
            Baseline = baseline;
            Baseline{name} = xvariables{name};
            S.(name) = simul_backward_model(initialconditions, lastperiod-firstperiod+1, Baseline);
        end
        % Compute marginal contributions (removing baseline)
        for j=1:length(variables)
            cumulatedcontribs = S.baseline{variables{j}}(firstperiod:lastperiod).data;
            contributions.(variables{j}) = dseries(transform(cumulatedcontribs), firstperiod, 'baseline');
            for i=1:xvariables.vobs
                name = xvariables.name{i};
                ts = S.(name);
                data = ts{variables{j}}(firstperiod:lastperiod).data;
                contributions.(variables{j}) = [contributions.(variables{j}), dseries(transform(data)-transform(cumulatedcontribs), firstperiod, name)];
            end
            contributions.(variables{j}) = contributions.(variables{j})(firstperiod:lastperiod);
        end
      otherwise
        error('Unknown method (%s)', method)
    end

    % Save output in caller workspace
    oname = getoutputname(varargin);
    assignin('caller', oname, contributions)

    % Cleanup
    rmdir(modfilename(1:end-4), 's')
    rmdir(sprintf('+%s', modfilename(1:end-4)), 's')
    delete(sprintf('%s.mod', modfilename(1:end-4)))
    delete(sprintf('%s.log', modfilename(1:end-4)))
end


function model = getmodel(cellarray)

    % Return variables for which we want to compute the contributions.
    %
    % INPUTS
    % - cellarray     [char]      1×n cell array of row char arrays.
    %
    % OUTPUTS
    % - var           [char]      name of the model (with extension)

    mpos = positions(cellarray);

    model = cellarray{mpos+1};

end


function eqtags = geteqtags(cellarray)

% Return equation tags for the equations we want to compute the contributions.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - eqtags        [char]      1×p cell array of row char arrays.

    [~, vpos, ~, ~, ~, ~, ~, ~, indices] = positions(cellarray);

    lastvalue = indices(find(indices==vpos)+1)-1;

    eqtags = cellarray(vpos+1:lastvalue);

end


function dname = getdatasetname(cellarray)

% Return the name of the dataset.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - dname         [char]      dataset name for endogenous and exogenous variables

    [~, ~, dpos] = positions(cellarray);

    dname = cellarray{dpos+1};

end


function [firstperiod, lastperiod] = getperiods(cellarray)

% Return variables for which we want to compute the contributions.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - ds            [dseries]   dataset for endogenous and exogenous variables

    [~, ~, ~, rpos] = positions(cellarray);

    firstperiod = dates();
    lastperiod = dates();

    if ~isempty(rpos)
        try
            tmp = strsplit(cellarray{rpos+1},':');
            firstperiod = dates(tmp{1});
            lastperiod = dates(tmp{2});
        catch
            error('dcontrib:getperiods: Cannot convert the --range argument to dates objects.')
        end
        if lastperiod<=firstperiod
            error('dcontrib:getperiods: In --range A:B we must have B>A.')
        end
    end

end


function dname = getbaselinename(cellarray)

% Return the name of the dataset.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - dname         [char]      baseline name for endogenous and exogenous variables

    [~, ~, ~, ~, bpos] = positions(cellarray);

    dname = cellarray{bpos+1};

end


function oname = getoutputname(cellarray)

% Return the name of the output.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - dname         [char]      baseline name for endogenous and exogenous variables

    [~, ~, ~, ~, ~, opos] = positions(cellarray);

    oname = cellarray{opos+1};

end


function method = getmethod(cellarray)

% Return the method for computing the dynaamic contributions.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - method        [char]      method: 'cumulate' or 'diff'

    [~, ~, ~, ~, ~, ~, kpos] = positions(cellarray);

    if isempty(kpos)
        method = 'cumulate';
    else
        method = cellarray{kpos+1};
    end

end


function bool = islog(cellarray)

% Returns true if the contributions are required in logs.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - method        [char]      method: 'cumulate' or 'diff'

    [~, ~, ~, ~, ~, ~, ~, lpos] = positions(cellarray);

   bool = ~isempty(lpos);

end



function [mpos, vpos, dpos, rpos, bpos, opos, kpos, lpos, indices] = positions(cellarray)

% Return  positions of the arguments.
%
% INPUTS
% - cellarray     [char]      1×n cell array of row char arrays.
%
% OUTPUTS
% - mpos          [integer]   scalar, index for the --model argument.
% - vpos          [integer]   scalar, index for the --tags arguments.
% - dpos          [integer]   scalar, index for the --database argument.
% - rpos          [integer]   scalar, index for the --range argument.
% - bpos          [integer]   scalar. index for the --baseline argument.
% - opos          [integer]   scalar, index for the --output argument.
% - kpos          [integer]   scalar, index for the --method argument.
% - lpos          [integer]   scalar, index for the --log option.

    % Index for --model argument
     mpos = find(strcmp('--model', cellarray));
     if isempty(mpos)
        error('dcontrib::positions: --model argument is mandatory.')
    elseif length(mpos)>1
        error('dplot::positions: Only one --model argument is allowed.')
    end

    % Index for --tags argument
    vpos = find(strcmp('--tags', cellarray));
    if isempty(vpos)
        error('dplot::positions: --tags argument is mandatory.')
    elseif length(vpos)>1
        error('dplot::positions: Only one --tags argument is allowed.')
    end

    % Index for the --initialconditions argument
     dpos = find(strcmp('--database', cellarray));
     if isempty(dpos)
        error('dplot::positions: --database argument is mandatory.')
    elseif length(dpos)>1
        error('dplot::positions: Only one --database argument is allowed.')
    end

    % Index for the --range argument
    rpos = find(strcmp('--range', cellarray));
    if length(rpos)>1
        error('dplot::positions: Only one --range argument is allowed.')
    end

    % Index for the --baseline argument
    bpos = find(strcmp('--baseline', cellarray));
    if length(bpos)>1
        error('dplot::positions: Only one --baseline argument is allowed.')
    end

    if ~isempty(rpos) && ~isempty(bpos)
        error('dplot::positions: --baseline and --range arguments are not allowed simultaneously.')
    end

    % Index for the --output argument.
    opos = find(strcmp('--output', cellarray));
    if isempty(opos)
        error('dplot::positions: --output argument is mandatory.')
    elseif length(opos)>1
        error('dplot::positions: Only one --periods argument is allowed.')
    end

    % Index for --method argument
    kpos = find(strcmp('--method', cellarray));
    if length(kpos)>1
        error('dplot::positions: Only one --method argument is allowed.')
    end

    % Index for --log option
     lpos = find(strcmp('--log', cellarray));
    if length(lpos)>1
        warning('dplot::positions: There is no point in using --log more than once.')
    end

    % Sorted vector of indices
     indices = sort([mpos; vpos; dpos; rpos; bpos; opos; kpos; lpos]);

end
