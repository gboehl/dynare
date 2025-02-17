function delete_dir_if_exists(dirname)
% function delete_dir_if_exists(dirname)
% removes directory if it exists
%
% INPUTS
%    dirname: directory to delete
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2011-2023 Dynare Team
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

if isoctave
    confirm_recursive_rmdir(false, 'local')
end

if exist(dirname, 'dir')
    rmdir(dirname, 's');
end
end
