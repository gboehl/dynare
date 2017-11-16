function [deviations, baseline, irfs] = backward_model_irf(initialcondition, innovationbaseline, listofshocks, listofvariables, varargin)

% Returns impulse response functions.
%
% INPUTS
% - initialcondition    [dseries]                      Initial conditions for the endogenous variables, or period 0.
% - innovationbaseline  [dseries]                      Baseline for the future innovations. If empty the baseline scenario is zero for future shocks.
% - listofshocks        [cell of strings or dseries]   The innovations for which the IRFs need to be computed.
% - listofvariables     [cell of strings]              The endogenous variables which will be returned.
% - periods             [integer]                      scalar, the number of periods.
%
% OUTPUTS
% - irfs                [struct of dseries]
%
% REMARKS
% - The names of the fields in the returned structure are given by the name
%   of the innovations listed in the second input argument. Each field gather
%   the associated paths for endogenous variables listed in the third input
%   argument.
% - If second argument is not empty, periods must not be greater than innovationbaseline.nobs.

% Copyright (C) 2017 Dynare Team
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

global M_ options_ oo_

% Check that the model is actually backward
if M_.maximum_lead
    error(['simul_model_irf:: The specified model is not backward looking!'])
end

% Set default value for the fourth input argument.
if nargin<5
    periods = 40;
    notransform = true;
else
    periods = varargin{1};
end

% Set default value for the last input argument (no transformation).
if nargin<6
    notransform = true;
else
    notransform = false;
    transform = varargin{2};
end

% Check third argument.
if ~iscell(listofshocks)
    error('Third input argument has to be a cell of string or dseries objects!')
else
    if all(cellfun(@ischar, listofshocks))
        deterministicshockflag = false;
    elseif all(cellfun(@isdseries, listofshocks))
        deterministicshockflag = true;
    else
        error('Elements of third input argument must all be char arrays or dseries objects!')
    end
    if deterministicshockflag
        numberofexperiments = length(listofshocks);
        exonames = cellstr(M_.exo_names);
        initialconditionperiod = initialcondition.dates(end);
        for i=1:numberofexperiments
            shock = listofshocks{i};
            impulsenames = shock.name;
            listofunknownexo = setdiff(impulsenames, exonames);
            if ~isempty(listofunknownexo)
                disp(listofunknownexo)
                error('In experiment n°%s, some of the declared shocks are unknown!', int2str(i))
            end
            if initialconditionperiod>=shock.dates(1)
                error('In experiment n°%s, the shock period must follow %s!', string(initialconditionperiod))
            end
        end
    end
end

baselineflag = false;
% Set default values for the baseline paths.
%
% TODO zero for all variables is probably a poor choice. It should be
% zero for additive exogenous variables and 1 for multiplicative
% exogenous variables.
Innovations = zeros(periods, M_.exo_nbr);

if ~isempty(innovationbaseline)
    if ~isdseries(innovationbaseline)
        error('If not empty, the second argument has to be a dseries object!')
    end
    if ~isequal(innovationbaseline.dates(1)-initialcondition.dates(end), 1)
        error('The first date of the second input argument must follow the last date of the first input argument!')
    end
    if innovationbaseline.nobs<periods
        error('The second input argument must at least have %s observations or lower the number of periods.', periods)
    end
    % Fill innovations with provided paths for the innovations.
    exonames = cellstr(M_.exo_names);
    for i = 1:length(exonames)
        if ~isempty(strmatch(exonames{i}, innovationbaseline.name))
            Innovations(:,i) = innovationbaseline{exonames{i}}.data(1:periods);
        end
    end
    baselineflag = true;
end

% Set up initial conditions
[initialcondition, periods, Innovations, DynareOptions, DynareModel, DynareOutput, endonames, exonames, nx, ny1, iy1, jdx, model_dynamic, y] = ...
    simul_backward_model_init(initialcondition, periods, options_, M_, oo_, Innovations);

% Get the covariance matrix of the shocks.
if ~deterministicshockflag
    if nnz(M_.Sigma_e)
        Sigma = M_.Sigma_e + 1e-14*eye(M_.exo_nbr);
        sigma = transpose(chol(Sigma));
    else
        error('You did not specify the size of the shocks!')
    end
end

% Initialization of the returned argument. Each will be a dseries object containing the IRFS for the endogenous variables listed in the third input argument.
deviations = struct();
baseline = dseries();
irfs = struct();

% Baseline paths (get transition paths induced by the initial condition and
% baseline innovations).
if options_.linear
    ysim__0 = simul_backward_linear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, Innovations, nx, ny1, iy1, jdx, model_dynamic);
else
    ysim__0 = simul_backward_nonlinear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, Innovations, iy1, model_dynamic);
end

% Transform the endogenous variables.
if notransform
    endo_simul__0 = ysim__0;
else
    endo_simul__0 = feval(transform, ysim__0);
end

% Compute the IRFs (loop over innovations).
for i=1:length(listofshocks)
    innovations = Innovations;
    % Add the shock.
    if deterministicshockflag
        shock = listofshocks{i};
        timid = shock.dates-initialconditionperiod;
        for j=1:shock.vobs
            k = find(strcmp(shock.name{i}, exonames));
            for l=1:length(timid)
                innovations(timid(l),k) = innovations(timid(l),k) + shock.data(l,j);
            end
        end
    else
        j = find(strcmp(listofshocks{i}, exonames));
        if isempty(j)
            error('backward_model_irf: Exogenous variable %s is unknown!', listofshocks{i})
        end
        innovations(1,:) = innovations(1,:) + transpose(sigma(:,j));
    end
    if options_.linear
        ysim__1 = simul_backward_linear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, nx, ny1, iy1, jdx, model_dynamic);
    else
        ysim__1 = simul_backward_nonlinear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, iy1, model_dynamic);
    end
    % Transform the endogenous variables
    if notransform
        endo_simul__1 = ysim__1;
    else
        endo_simul__1 = feval(transform, ysim__1);
    end
    % Instantiate a dseries object (with all the endogenous variables)
    alldeviations = dseries(transpose(endo_simul__1-endo_simul__0), initialcondition.init, endonames(1:M_.orig_endo_nbr), cellstr(DynareModel.endo_names_tex(1:M_.orig_endo_nbr,:)));
    if nargout>2
        allirfs = dseries(transpose(endo_simul__1), initialcondition.init, endonames(1:M_.orig_endo_nbr), cellstr(DynareModel.endo_names_tex(1:M_.orig_endo_nbr,:)));
    end
    % Extract a sub-dseries object
    if deterministicshockflag
        name = sprintf('experiment_%s', int2str(i));
    else
        name = listofshocks{i};
    end
    deviations.(name) = alldeviations{listofvariables{:}};
    if nargout>2
        irfs.(name) = allirfs{listofvariables{:}};
        irfs.(name) = [irfs.(name) dseries(innovations, initialcondition.last+1, exonames)];
    end
end

if nargout>1
    baseline = dseries(transpose(endo_simul__0), initialcondition.init, endonames(1:M_.orig_endo_nbr), cellstr(DynareModel.endo_names_tex(1:M_.orig_endo_nbr,:)));
    baseline = [baseline, innovationbaseline];
end