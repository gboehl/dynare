function dummy = set(pacmodel, paceq, parameter, isnonzeromean)

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

global M_

eqtag = M_.pac.(pacmodel).tag_map{strcmp(paceq, M_.pac.(pacmodel).tag_map(:,1)),2};

dummy = [];
idp = find(strcmp(parameter, M_.param_names));
if ~isempty(idp)
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'additive') && length(M_.pac.(pacmodel).equations.(eqtag).additive.vars)>1
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).additive, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).additive.bgp = false(1, length(M_.pac.(pacmodel).equations.(eqtag).additive.params));
        end
        [isparam, ip] = ismember(idp, M_.pac.(pacmodel).equations.(eqtag).additive.params);
        if isparam
            M_.pac.(pacmodel).equations.(eqtag).additive.bgp(ip) = isnonzeromean;
            return
        end
    end
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'optim_additive')
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).optim_additive, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).optim_additive.bgp = false(1, length(M_.pac.(pacmodel).equations.(eqtag).optim_additive.params));
        end
        [isparam, ip] = ismember(idp, M_.pac.(pacmodel).equations.(eqtag).optim_additive.params);
        if isparam
            M_.pac.(pacmodel).equations.(eqtag).optim_additive.bgp(ip) = isnonzeromean;
            return
        end
    end
    if isfield(M_.pac.(pacmodel).equations.(eqtag), 'non_optimizing_behaviour')
        if ~isfield(M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour, 'bgp')
            M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.bgp = false(1, length(M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.params));
        end
        [isparam, ip] = ismember(idp, M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.params);
        if isparam
            M_.pac.(pacmodel).equations.(eqtag).non_optimizing_behaviour.bgp(ip) = isnonzeromean;
            return
        end
    end
    warning('Parameter %s is not associated to an exogenous variable in equation %s.', parameter, paceq)
else
    error('Parameter %s is unknown.', parameter)
end