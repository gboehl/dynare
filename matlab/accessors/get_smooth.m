
function y0 = get_smooth(varargin)
% function x = get_smooth(vname1, vname2, )
% returns smoothed variables or shocks identified by their name
%
% INPUTS:
%   vname1, vname2, ... :  list of variable/shock names
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

SmoothedVariables=[struct2cell(oo_.SmoothedVariables); struct2cell(oo_.SmoothedShocks)];
my_field_names = [fieldnames(oo_.SmoothedVariables); fieldnames(oo_.SmoothedShocks)];
isvar=zeros(length(SmoothedVariables),1);
for jf = 1:length(SmoothedVariables)
    isvar(jf)=~(isstruct(SmoothedVariables{jf}));
end
SmoothedVariables=cell2struct(SmoothedVariables(logical(isvar)),my_field_names(logical(isvar)));


y0=zeros(length(SmoothedVariables.(varargin{1})),length(varargin));
for j=1:length(varargin)
    y0(:,j)=SmoothedVariables.(varargin{j});
end
