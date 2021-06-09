function [Yvec, lhssubvec, Xmat, constrained] = put_in_sur_form(Y, lhssub, X)
%function [Yvec, lhssubvec, Xmat, constrained] = put_in_sur_form(Y, lhssub, X)
%
% INPUTS
%   Y           [cell array]  dependent variables
%   lhssub      [cell array]  RHS subtracted from LHS
%   X           [cell array]  regressors
%
% OUTPUTS
%   Yvec        [vector]      dependent variables
%   lhssubvec   [cell array]  RHS subtracted from LHS
%   Xmat        [matrix]      regressors
%   constrained [cellstr]     names of parameters that were constrained
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

%% Check inputs
if nargin ~= 3
    error('function expects 3 arguments');
end

if isempty(Y) || ~iscell(Y) ...
        || isempty(lhssub) || ~iscell(lhssub) ...
        || isempty(X) || ~iscell(X) ...
        || length(Y) ~= length(X) ...
        || length(Y) ~= length(lhssub)
    error('function arguments should be cells of the same size');
end

%% Organize output
neqs = length(Y);
nobs = zeros(neqs, 1);
for i = 1:neqs
    if ~isempty(X{i})
        % X{i} is empty for AR(1) equations
        assert(size(X{i}, 1) == size(Y{i}, 1), 'Y{i} and X{i} must have the same nuber of observations');
    end
    nobs(i) = size(Y{i}, 1);
end
fd = Y{1}.firstdate;
nrows = sum(nobs);
Xmat = dseries();
Yvec = dseries();
lhssubvec = dseries();
constrained = {};
for i = 1:neqs
    if ~isempty(X{i})
        to_remove = [];
        nr = sum(nobs(1:i-1));
        nxcol = size(X{i}, 2);
        names = X{i}.name;
        Xtmp = dseries([zeros(nr, nxcol); X{i}.data; zeros(nrows-nr-nobs(i), nxcol)], fd, names);
        Xmatnames = Xmat.name;
        for j = 1:length(names)
            idx = find(strcmp(Xmatnames, names{j}));
            if ~isempty(idx)
                Xmat.(Xmatnames{idx}) = Xmat.(Xmatnames{idx}) + Xtmp.(names{j});
                Xtmp = Xtmp.remove(names{j});
                constrained{end+1} = Xmatnames{idx};
            end
        end
        if ~isempty(Xtmp)
            Xmat = [Xmat Xtmp];
        end
    end
    Yvec = dseries([Yvec.data; Y{i}.data], fd);
    if isempty(lhssub{i})
        lhssubvec = dseries([lhssubvec.data; zeros(size(Y{i}.data, 1), 1)], fd);
    else
        lhssubvec = dseries([lhssubvec.data; lhssub{i}.data], fd);
    end
end
assert(size(Yvec, 1) == size(Xmat, 1));
assert(size(Yvec, 1) == size(lhssubvec, 1));
end
