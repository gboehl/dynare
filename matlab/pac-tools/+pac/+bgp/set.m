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
if ~isempty(ide)
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'additive') && length(M_.pac.(pacmodel).equations.(eqtag).additive.vars)>1
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).additive, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).additive.bgp = zeros(1, length(M_.pac.(pacmodel).equations.(eqtag).additive.params));
        end
        [isvar, ie] = ismember(ide, M_.pac.(pacmodel).equations.(eqtag).additive.vars);
        if isvar
            M_.pac.(pacmodel).equations.(eqtag).additive.bgp(ie) = nonzeromean;
            return
        end
    end
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'optim_additive')
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).optim_additive, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).optim_additive.bgp = zeros(1, length(M_.pac.(pacmodel).equations.(eqtag).optim_additive.params));
        end
        [isvar, ie] = ismember(ide, M_.pac.(pacmodel).equations.(eqtag).optim_additive.vars);
        if isvar
            M_.pac.(pacmodel).equations.(eqtag).optim_additive.bgp(ie) = nonzeromean;
            return
        end
    end
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'non_optimizing_behaviour')
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.bgp = zeros(1, length(M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.params));
        end
        [isvar, ie] = ismember(ide, M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.vars);
        if isvar
            M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.bgp(ie) = nonzeromean;
            return
        end
    end
    warning('%s is not an exogenous variable in equation %s.', variable, paceq)
else
    error('Endogenous variable %s is unknown.', variable)
end