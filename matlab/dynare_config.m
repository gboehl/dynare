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

% Copyright © 2001-2023 Dynare Team
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

if nargin && ~isempty(path_to_dynare)
    addpath(path_to_dynare);
end

dynareroot = strrep(which('dynare'),'dynare.m','');

origin = pwd();
cd([dynareroot '/..'])

p = {'/../contrib/ms-sbvar/TZcode/MatlabFiles/' ; ...
     '/accessors/' ; ...
     '/AIM/' ; ...
     '/backward/' ; ...
     '/cli/' ; ...
     '/convergence_diagnostics/' ; ...
     '/discretionary_policy/' ; ...
     '/distributions/' ; ...
     '/ep/' ; ...
     '/estimation/'; ...
     '/estimation/smc/'; ...
     '/estimation/resampler/'; ...
     '/kalman/' ; ...
     '/kalman/likelihood' ; ...
     '/latex/' ; ...
     '/lmmcp/' ; ...
     '/dseries/src/' ; ...
     '/reporting/' ; ...
     '/matrix_solver/'; ...
     '/moments/'; ...
     '/ms-sbvar/' ; ...
     '/ms-sbvar/identification/' ; ...
     '/nonlinear-filters/' ; ...
     '/ols/' ; ...
     '/optimal_policy/' ; ...
     '/optimization/' ; ...
     '/pac-tools/' ; ...
     '/parallel/' ; ...
     '/partial_information/' ; ...
     '/perfect-foresight-models/' ; ...
     '/shock_decomposition/' ; ...
     '/stochastic_solver/' ; ...
     '/utilities/dataset/' ; ...
     '/utilities/doc/' ; ...
     '/utilities/estimation/' ; ...
     '/utilities/general/' ; ...
     '/utilities/graphics/' ; ...
     '/utilities/tests/src/' ; ...
     '/utilities/version/'};

% For functions that exist only under some Octave versions
% or some MATLAB versions, and for which we provide some replacement functions

% Replacements for rows(), columns(), vec() and issquare() (inexistent under MATLAB)
if ~isoctave
    p{end+1} = '/missing/rows_columns';
    p{end+1} = '/missing/issquare';
    p{end+1} = '/missing/vec';
end

% Replacements for functions of the MATLAB statistics toolbox
if isoctave
    % Under Octave, these functions are in the statistics Forge package.
    % Our replacement functions don't work under Octave (because of gamrnd, see
    % #1638), hence the statistics toolbox is now a hard requirement
    if ~user_has_octave_forge_package('statistics')
        error('You must install the "statistics" package from Octave Forge, either with your distribution package manager or with "pkg install -forge io statistics"')
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

% contains and splitlines don't exist in Octave
if isoctave
    p{end+1} = '/missing/contains';
    p{end+1} = '/missing/splitlines';
end

% datetime doesn't exist in Octave
if isoctave
    p{end+1} = '/missing/datetime';
end

% intersect with 'stable' flag is broken before Octave 8.4, bug #60347
if isoctave && octave_ver_less_than('8.4')
    p{end+1} = '/missing/intersect_stable';
end

P = cellfun(@(c)[dynareroot(1:end-1) c], p, 'uni',false);

% Get mex files folder(s)
mexpaths = get_path_to_mex_files(dynareroot);

% Add mex files folder(s)
P(end+1:end+length(mexpaths)) = mexpaths;

% Set matlab's path
addpath(P{:});

% Initialization of the dates and dseries classes (recursive).
initialize_dseries_class();

cd(origin);
