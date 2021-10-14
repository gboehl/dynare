function dummy = set(pacmodel, paceq, variable, nonzeromean)

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

eqtag = M_.pac.(pacmodel).tag_map{strcmp(paceq, M_.pac.(pacmodel).tag_map(:,1)),2};

dummy = [];

ide = find(strcmp(variable, M_.endo_names));
xflag = false;
if isempty(ide)
    % variable is not an endogenous variable
    ide = find(strcmp(variable, M_.exo_names));
    xflag = true;
end

if ~isempty(ide)
    [M_, done] = writebgpfield('additive', pacmodel, eqtag, ide, xflag, nonzeromean, M_);
    if done, return, end
    [M_, done] = writebgpfield('optim_additive', pacmodel, eqtag, ide, xflag, nonzeromean, M_);
    if done, return, end
    [M_, done] = writebgpfield('non_optimizing_behaviour', pacmodel, eqtag, ide, xflag, nonzeromean, M_);
    if done, return, end
    warning('%s is not an exogenous variable in equation %s.', variable, paceq)
else
    error('Endogenous/Exogenous variable %s is unknown.', variable)
end


function [M_, done] = writebgpfield(type, pacmodel, eqtag, ide, xflag, nonzeromean, M_)
done = false;
if isfield(M_.pac.(pacmodel).equations.(eqtag), type)
    if ~isfield(M_.pac.(pacmodel).equations.(eqtag).additive, 'bgp')
        M_.pac.(pacmodel).equations.(eqtag).(type).bgp = repmat({false}, 1, length(M_.pac.(pacmodel).equations.(eqtag).(type).params));
    end
    [isvar, ie] = ismember(ide, M_.pac.(pacmodel).equations.(eqtag).(type).vars);
    if isvar
        if xflag
            assert(~M_.pac.(pacmodel).equations.(eqtag).(type).isendo(ie), 'Variable type issue.')
        else
            assert(M_.pac.(pacmodel).equations.(eqtag).(type).isendo(ie), 'Variable type issue.')
        end
        M_.pac.(pacmodel).equations.(eqtag).(type).bgp{ie} = nonzeromean;
        done = true;
    end
end
