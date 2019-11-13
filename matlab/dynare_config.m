function dynareroot = dynare_config(path_to_dynare)
%function dynareroot = dynare_config(path_to_dynare)
%
% This function tests the existence of valid mex files (for qz
% decomposition, solution to sylvester equation and kronecker
% products...) and, if needed, add paths to the matlab versions
% of these routines.
% Also adds other directories to the path.
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2019 Dynare Team
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

if nargin && ~isempty(path_to_dynare)
    addpath(path_to_dynare);
end

dynareroot = strrep(which('dynare'),'dynare.m','');

origin = pwd();
cd([dynareroot '/..'])

p = {'/distributions/' ; ...
     '/kalman/' ; ...
     '/kalman/likelihood' ; ...
     '/AIM/' ; ...
     '/partial_information/' ; ...
     '/perfect-foresight-models/' ; ...
     '/ms-sbvar/' ; ...
     '/ms-sbvar/identification/' ; ...
     '/../contrib/ms-sbvar/TZcode/MatlabFiles/' ; ...
     '/../contrib/jsonlab/' ; ...
     '/parallel/' ; ...
     '/particles/src' ; ...
     '/gsa/' ; ...
     '/ep/' ; ...
     '/backward/' ; ...
     '/convergence_diagnostics/' ; ... 
     '/cli/' ; ...
     '/lmmcp/' ; ...
     '/optimization/' ; ...
     '/modules/dseries/src/' ; ...
     '/utilities/doc/' ; ...
     '/utilities/tests/src/' ; ...
     '/utilities/dataset/' ; ...
     '/utilities/general/' ; ...
     '/utilities/graphics/' ; ...
     '/modules/reporting/src/' ; ...
     '/modules/reporting/macros/'};

% For functions that exist only under some Octave versions
% or some MATLAB versions, and for which we provide some replacement functions

% Replacements for rows(), columns(), vec() and issquare() (inexistent under MATLAB)
if ~isoctave
    p{end+1} = '/missing/rows_columns';
    p{end+1} = '/missing/issquare';
    p{end+1} = '/missing/vec';
end

% ordeig() doesn't exist in Octave < 5
if isoctave && octave_ver_less_than('5')
    p{end+1} = '/missing/ordeig';
end

% corrcoef with two outputs is missing in Octave < 4.4 (ticket #796)
if isoctave && octave_ver_less_than('4.4') && ~user_has_octave_forge_package('nan')
    p{end+1} = '/missing/corrcoef';
end

%% intersect(…, 'stable') doesn't exist in Octave and in MATLAB < R2013a
if isoctave || matlab_ver_less_than('8.1')
    p{end+1} = '/missing/intersect_stable';
end

% Replacements for functions of the MATLAB statistics toolbox
if isoctave
    % These functions were part of Octave < 4.4, they are now in the statistics Forge package
    if ~octave_ver_less_than('4.4') && ~user_has_octave_forge_package('statistics')
        % Our replacement functions don't work under Octave (because of gamrnd, see
        % #1638), hence the statistics toolbox is now a hard requirement
        error('You must install the "statistics" package from Octave Forge, either with your distribution package manager or with "pkg install -forge statistics"')
    end
else
    if ~user_has_matlab_license('statistics_toolbox')
        p{end+1} = '/missing/stats/';
    end
end

% Check if struct2array is available.
if ~exist('struct2array')
    p{end+1} = '/missing/struct2array';
end

% isfile is missing in Octave < 5 and Matlab < R2017b
if (isoctave && octave_ver_less_than('5')) || (~isoctave && matlab_ver_less_than('9.3'))
    p{end+1} = '/missing/isfile';
end

% strsplit and strjoin are missing in MATLAB < R2013a
if ~isoctave && matlab_ver_less_than('8.1')
    p{end+1} = '/missing/strsplit';
    p{end+1} = '/missing/strjoin';
end

% contains and splitlines don't exist in Octave and in MATLAB < R2016b
if isoctave || matlab_ver_less_than('9.1')
    p{end+1} = '/missing/contains';
    p{end+1} = '/missing/splitlines';
end

% isrow, iscolumn and ismatrix are missing in Matlab<R2010b
if ~isoctave && matlab_ver_less_than('7.11')
    p{end+1} = '/missing/is-row-column-matrix';
end

%% isdiag is missing in MATLAB < R2014a
if ~isoctave && matlab_ver_less_than('8.4')
    p{end+1} = '/missing/isdiag';
end

%% narginchk is missing in MATLAB < R2011b
if ~isoctave && matlab_ver_less_than('7.13')
    p{end+1} = '/missing/narginchk';
end

P = cellfun(@(c)[dynareroot(1:end-1) c], p, 'uni',false);

% Get mex files folder(s)
mexpaths = add_path_to_mex_files(dynareroot, false);

% Add mex files folder(s)
P(end+1:end+length(mexpaths)) = mexpaths;

% Set matlab's path
addpath(P{:});

% Initialization of the dates and dseries classes (recursive).
initialize_dseries_class();

cd(origin);
