function write_param_init_inc_file(subfolder, fnameroot, idxs, estimated_params)
%function write_param_init_inc_file(fname, idxs, estimated_params)
% Write .inc files that initialize parameters
%
% INPUTS
%   subfolder        [string] subfolder in which to place file
%   fnameroot        [string] root of filename
%   idxs             [vector] indexes in M_.params of estimated parameters
%   estimated_params [vector] estimated parameters  
%
% OUTPUTS
%   none
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

global M_

%% Check inputs
if nargin ~= 4
    error('function takes 4 arguments')
end

if ~ischar(subfolder)
    error([subfolder ' must be a string']);
end

if ~ischar(fnameroot)
    error([fnameroot ' must be a string']);
end

if ~isvector(idxs) ...
        || ~isvector(estimated_params) ...
        || length(idxs) ~= length(estimated_params)
    error('the two arguments must be vectors of the same length')
end

%% Write file
% Open
filepath = [M_.fname filesep 'model' filesep subfolder];
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

fname = [filepath filesep fnameroot '-declare.inc'];
fid = fopen(fname, 'w');
if fid < 0
    error(['couldn''t open file' fname])
end

% Write declaration
fprintf(fid, 'parameters');
for i = 1:length(idxs)
    fprintf(fid, ' %s', M_.param_names{idxs(i)});
end
fprintf(fid, ';\n');
fclose(fid);

fname = [filepath filesep fnameroot '-initialize.inc'];
fid = fopen(fname, 'w');
if fid < 0
    error(['couldn''t open file' fname])
end
for i = 1:length(idxs)
    fprintf(fid, '%s = %f;\n', M_.param_names{idxs(i)}, estimated_params(i));
end
fclose(fid);
end
