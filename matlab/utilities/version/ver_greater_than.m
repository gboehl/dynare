function tf = ver_greater_than(ver1, ver2)  % --*-- Unitary tests --*--
%function tf = ver_greater_than(ver1, ver2)
% ver1 > ver2 ? 1 : 0;
%
% INPUTS
%    ver1    [string]    software version number
%    ver2    [string]    software version number
%
% OUTPUTS
%    tf      [bool]      true if ver1 > ver2
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2015-2021 Dynare Team
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

tf = ~ver_less_than(ver1, ver2) && ~strcmp(ver1, ver2);
return

%@test:1
ver2='4.4';
ver1='4.5.2';
t(1)=dassert(ver_greater_than(ver1,ver2),true);
T = all(t);
%@eof:1
%@test:2
ver1='6-unstable-2021-12-15-1737-21a8a579';
ver2='4.4';
t(1)=dassert(ver_greater_than(ver1,ver2),true);
T = all(t);
%@eof:2
%@test:3
ver2='5.0';
ver1='5.1';
t(1)=dassert(ver_greater_than(ver1,ver2),true);
T = all(t);
%@eof:3
%@test:4
ver2='6-unstable-2021-12-18-1227-c43777f6';
ver1='6-unstable-2021-12-19-1953-d841fc7c';
t(1)=dassert(ver_greater_than(ver1,ver2),true);
T = all(t);
%@eof:4