function B = olseqs(ds, varargin)
%function B = olseqs(ds, varargin)
% Run OLS on chosen model equations
%
% INPUTS
%   ds      [dseries]    data
%
% OUTPUTS
%   B       [vector]     estimated parameters
%
% SPECIAL REQUIREMENTS
%   none

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

global M_

assert(nargin == 1 || nargin == 3, 'Incorrect number of arguments passed to olseqs');

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, ~] = getEquationsByTags(jsonmodel, varargin{:});

Y = ds{lhs}.data;

rhs_ = strsplit(rhs, {'+','-','*','/','^','log(','exp(','(',')'});
rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
vnames = setdiff(rhs_, cellstr(M_.param_names));
regexpr = cell2mat(strcat('(', vnames, {'\(-\d+\))|'}));
vwlags = regexp(rhs, regexpr(1:end-1), 'match');
X = cell2mat(cellfun(@eval, strcat('ds.', vwlags, '.data'), 'UniformOutput', false));

% Remove all rows that have a NaN
[row, ~] = find(isnan(X), 1, 'last');
Y = Y(row+1:end, :);
X = X(row+1:end, :);

% OLS Estimation
B = (X'*X)\X'*Y;

end