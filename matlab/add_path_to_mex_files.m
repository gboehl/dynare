function mexpath = add_path_to_mex_files(dynareroot, modifypath)

% Copyright (C) 2015-2022 Dynare Team
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

if nargin<2
    modifypath = true;
end

if isoctave
    % Add specific paths for Dynare Windows package
    if ispc
        if strcmpi(computer(), 'i686-w64-mingw32')
            warning('MEX files not available for 32-bit Octave')
        else
            tmp = [dynareroot '../mex/octave/win64/'];
            if exist(tmp, 'dir')
                mexpath = tmp;
                if modifypath
                    addpath(mexpath);
                end
            end
        end
    end
    % Add generic MATLAB path (with higher priority than the previous ones)
    if exist('mexpath')
        mexpath = { mexpath; [dynareroot '../mex/octave/'] };
    else
        mexpath = { [dynareroot '../mex/octave/'] };
    end
    if modifypath
        addpath([dynareroot '../mex/octave/']);
    end
else
    if strcmp(computer, 'PCWIN')
        warning('MEX files not available for 32-bit MATLAB')
    end
    % Add win64 specific paths for Dynare Windows package
    if strcmp(computer, 'PCWIN64')
        if matlab_ver_less_than('9.4')
            tmp = [dynareroot '../mex/matlab/win64-8.3-9.3/'];
            if exist(tmp, 'dir')
                mexpath = tmp;
                if modifypath
                    addpath(mexpath);
                end
            end
        else
            tmp = [dynareroot '../mex/matlab/win64-9.4-9.12/'];
            if exist(tmp, 'dir')
                mexpath = tmp;
                if modifypath
                    addpath(mexpath);
                end
            end
        end
    end
    % Add macOS paths for Dynare Mac package
    if strcmp(computer, 'MACI64')
        if matlab_ver_less_than('9.4')
            tmp = [dynareroot '../mex/matlab/maci64-8.3-9.3/'];
            if exist(tmp, 'dir')
                mexpath = tmp;
                if modifypath
                    addpath(mexpath);
                end
            end
        else
            tmp = [dynareroot '../mex/matlab/maci64-9.4-9.12'];
            if exist(tmp, 'dir')
                mexpath = tmp;
                if modifypath
                    addpath(mexpath);
                end
            end
        end
    end
    % Add generic MATLAB path (with higher priority than the previous ones)
    if exist('mexpath')
        mexpath = { mexpath; [dynareroot '../mex/matlab/'] };
    else
        mexpath = { [dynareroot '../mex/matlab/'] };
    end
    if modifypath
        addpath([dynareroot '../mex/matlab/']);
    end
end
