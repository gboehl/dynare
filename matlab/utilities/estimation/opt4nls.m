function [noprint, opt] = opt4nls(varargin)

% Sets options for NLS routines.

% Copyright © 021 Dynare Team
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

varargin = varargin{1};
nargin = length(varargin);

noprint = false;

if nargin
    if mod(nargin, 2)
        error('Options must come by key/value pairs.')
    end
    i = 1;
    while i<nargin
        if isequal(varargin{i}, 'noprint')
            noprint = varargin{i+1};
            i = i+2;
            continue
        else
            if ~exist('opt', 'var')
                opt = sprintf('''%s''', varargin{i});
            else
                opt = sprintf('%s,''%s''', opt, varargin{i});
            end
            if isnumeric(varargin{i+1})
                opt = sprintf('%s,%s', opt, num2str(varargin{i+1}));
            else
                opt = sprintf('%s,''%s''', opt, varargin{i+1});
            end
            i = i+2;
        end
    end
else
    opt = '';
end