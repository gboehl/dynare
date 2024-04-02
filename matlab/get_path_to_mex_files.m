function mexpath = get_path_to_mex_files(dynareroot)
% Returns a cell array containing one or several directory paths
% which should contain the MEX files.

% Copyright © 2015-2024 Dynare Team
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

build_dir = get_build_dir(dynareroot);
if ~isempty(build_dir)
    % If a Meson build directory is found, use it preferably
    mexpath = { build_dir };
elseif isoctave
    % Add specific paths for Dynare Windows package
    if ispc
        if strcmpi(computer(), 'i686-w64-mingw32')
            warning('MEX files not available for 32-bit Octave')
        else
            tmp = [dynareroot '../mex/octave/win64/'];
            if exist(tmp, 'dir')
                mexpath = tmp;
            end
        end
    end
    % Add generic Octave path (with higher priority than the previous ones)
    if exist('mexpath')
        mexpath = { mexpath; [dynareroot '../mex/octave/'] };
    else
        mexpath = { [dynareroot '../mex/octave/'] };
    end
else
    if strcmp(computer, 'PCWIN')
        warning('MEX files not available for 32-bit MATLAB')
    end
    % Add win64 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN64')
        tmp = [dynareroot '../mex/matlab/win64-9.5-24.1/'];
        if exist(tmp, 'dir')
            mexpath = tmp;
        end
    end
    % Add macOS paths for Dynare Mac package
    if strcmp(computer, 'MACI64')
        tmp = [dynareroot '../mex/matlab/maci64-9.5-24.1/'];
        if exist(tmp, 'dir')
            mexpath = tmp;
        end
    end
    if strcmp(computer, 'MACA64')
        tmp = [dynareroot '../mex/matlab/maca64-23.2-24.1/'];
        if exist(tmp, 'dir')
            mexpath = tmp;
        end
    end
    % Add generic MATLAB path (with higher priority than the previous ones)
    if exist('mexpath')
        mexpath = { mexpath; [dynareroot '../mex/matlab/'] };
    else
        mexpath = { [dynareroot '../mex/matlab/'] };
    end
end
