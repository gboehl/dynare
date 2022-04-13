function estimate(method, data, varargin)
    
% Copyright © 2017 Dynare Team
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

% Check first input.
if ~ischar(method) || isempty(regexp(method, 'method\(((?=[\w])[^_0-9][\w]*)\)', 'match'))
    error('estimate:: First argument must be an estimation method')
end

% Check second input argument.
if ~ischar(data) || isempty(regexp(data, 'data\(((?=[\w])[^_0-9][\w]*)\)', 'match'))
    error('estimate:: Second argument must be a data option.')
end

% Get the estimation method.
tmp = regexp(method, 'method\(((?=[\w])[^_0-9][\w]*)\)', 'tokens');
method = tmp{1}{1};

% Get the data.
tmp = regexp(data, 'data\(((?=[\w])[^_0-9][\w]*)\)', 'tokens');
ds = tmp{1}{1};
if ismember(ds, evalin('caller','who'))
    ts = evalin('caller', ds); 
    if ~isdseries(ts)
        error('estimate:: %s has to be a dseries object!', ds)
    end
else
    error('estimate:: %s is unknown!', ds)
end

eqns = varargin(:);
nqns = length(varargin);

% Run estimations
for i=1:nqns
    if ~ismember(eqns{i}, M_.equations_tags(strmatch('name', M_.equations_tags(:,2)),3))
        error('estimate:: There is no equation named as %s!', eqns{i})
    end
    switch method
      case 'ols'
        olseqs(ts, eqns{i});
      otherwise
        error('estimate:: Unknown estimation method')
    end
end