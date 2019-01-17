function [Yvec, Xmat, constrained] = put_in_sur_form(Y, X)
%function [Yvec, Xmat, constrained] = put_in_sur_form(Y, X)
%
% INPUTS
%   Y           [cell array]  dependent variables
%   X           [cell array]  regressors
%
% OUTPUTS
%   Yvec        [vector]      dependent variables
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%% Check inputs
if nargin ~= 2
    error('put_in_sur_form expects 2 arguments');
end

if isempty(Y) || ~iscell(Y) || isempty(X) || ~iscell(X) || length(Y) ~= length(X)
    error('put_in_sur_form arguments should be cells of the same size');
end

%% Organize output
neqs = length(Y);
nobs = zeros(neqs, 1);
for i = 1:neqs
    if ~isempty(X{i})
        % X{i} is empty for AR(1) equations
        assert(size(X{i}, 1) == size(Y{i}, 1), 'Y{i} and X{i} must have the same nuber of observations');
    end
    nobs(i) = size(X{i}, 1);
end
nrows = sum(nobs);
Xmat = dseries([X{1}.data; zeros(nrows-nobs(1), size(X{1}, 2))], X{1}.firstdate, X{1}.name);
Yvec = Y{1};
constrained = {};
for i = 2:neqs
    if isempty(X{i})
        data = [Xmat.data;zeros(size(Y{i},1), size(Xmat.data, 2))];
        Xmat = dseries(data, X{1}.firstdate, Xmat.name);
    else
        to_remove = [];
        nr = sum(nobs(1:i-1));
        nxcol = size(X{i}, 2);
        Xtmp = dseries([zeros(nr, nxcol); X{i}.data; zeros(nrows-nr-nobs(i), nxcol)], X{1}.firstdate, X{i}.name);
        for j = 1:length(X{i}.name)
            idx = find(strcmp(Xmat.name, X{i}.name{j}));
            if ~isempty(idx)
                Xmat.(Xmat.name{idx}) = Xmat{idx} + Xtmp{j};
                to_remove = [to_remove j];
                constrained{end+1} = Xmat.name{idx};
            end
        end
        for j = length(to_remove):-1:1
            Xtmp = Xtmp.remove(Xtmp.name{j});
        end
        if ~isempty(Xtmp)
            Xmat = [Xmat Xtmp];
        end
    end
    Yvec = dseries([Yvec.data; Y{i}.data], Yvec.firstdate);
end
end
