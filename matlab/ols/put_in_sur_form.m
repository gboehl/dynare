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
nobs = size(X{1}, 1);
neqs = length(X);
Xmat = dseries([X{1}.data; zeros(nobs*(neqs-1), size(X{1}, 2))], X{1}.firstdate, X{1}.name);
Yvec = Y{1};
constrained = {};
for i = 2:neqs
    to_remove = [];
    Xtmp = dseries([zeros(nobs*(i-1), size(X{i}, 2)); X{i}.data; zeros(nobs*(neqs-i), size(X{i}, 2))], X{i}.firstdate, X{i}.name);
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
    Yvec = dseries([Yvec.data; Y{i}.data], Yvec.firstdate);
end
end
