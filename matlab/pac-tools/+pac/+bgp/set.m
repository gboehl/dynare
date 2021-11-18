function dummy = set(pacmodl, paceq, variable, nonzeromean)

% Provide information about long run levels of exogenous variables in PAC equation.
%
% INPUTS
% - pacmodel       [char]    1×n array, name of the PAC model
% - paceq          [char]    1×m array, name of the PAC equation
% - variable       [char]    1×p array, name of the variable (exogenous variable in the PAC equation)
% - nonzeromean    [double]  scalar, mean of the exogenous variable
%
% OUPUTS
% - dummy          [double]  empty array.

% Copyright © 2019-2021 Dynare Team
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

dummy = [];

pacmodel = M_.pac.(pacmodl);

ide = find(strcmp(variable, M_.endo_names));
xflag = false;
if isempty(ide)
    % variable is not an endogenous variable
    ide = find(strcmp(variable, M_.exo_names));
    xflag = true;
end

if ~isempty(ide)
    [pacmodel, done] = writebgpfield('additive', pacmodel, ide, xflag, nonzeromean);
    if done, M_.pac.(pacmodl) = pacmodel; return, end
    [pacmodel, done] = writebgpfield('optim_additive', pacmodel, ide, xflag, nonzeromean);
    if done, M_.pac.(pacmodl) = pacmodel; return, end
    [pacmodel, done] = writebgpfield('non_optimizing_behaviour', pacmodel, ide, xflag, nonzeromean);
    if done, M_.pac.(pacmodl) = pacmodel; return, end
    warning('%s is not an exogenous variable in equation %s.', variable, paceq)
else
    error('Endogenous/Exogenous variable %s is unknown.', variable)
end


function [pacmodel, done] = writebgpfield(type, pacmodel, ide, xflag, nonzeromean, M_)
done = false;
if isfield(pacmodel, type)
    if ~isfield(pacmodel.additive, 'bgp')
        pacmodel.(type).bgp = repmat({false}, 1, length(pacmodel.(type).params));
    end
    [isvar, ie] = ismember(ide, pacmodel.(type).vars);
    if isvar
        if xflag
            assert(~pacmodel.(type).isendo(ie), 'Variable type issue.')
        else
            assert(pacmodel.(type).isendo(ie), 'Variable type issue.')
        end
        pacmodel.(type).bgp{ie} = nonzeromean;
        done = true;
    end
end