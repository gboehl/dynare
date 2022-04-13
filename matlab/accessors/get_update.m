function y0 = get_update(varargin)
% function x = get_update(vname1, vname2, )
% returns updated variables identified by their name
%
% INPUTS:
%   vname1, vname2, ... :  list of variable names
%
% OUTPUTS
%   x:      smoothed variables [T x number of variables]
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019 Dynare Team
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
global oo_

y0=zeros(length(oo_.UpdatedVariables.(varargin{1})),length(varargin));
for j=1:length(varargin)
    y0(:,j)=oo_.UpdatedVariables.(varargin{j});
end
