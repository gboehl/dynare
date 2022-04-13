function warning_config()
% Activates useful warnings
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2008-2022 Dynare Team
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

warning on;

% Display a calling stack trace when a warning is issued
warning('on', 'backtrace');

if isoctave
    warning('off', 'Octave:separator-insert');
    warning('off', 'Octave:language-extension');
    warning('off', 'Octave:single-quote-string');
    warning('off', 'Octave:missing-semicolon');
    warning('off', 'Octave:empty-list-elements');
    warning('off', 'Octave:num-to-str');
    warning('off', 'Octave:resize-on-range-error');
    warning('off', 'Octave:str-to-num');
    warning('off', 'Octave:array-to-scalar');
    warning('off', 'Octave:array-to-vector');
    warning('off', 'Octave:mixed-string-concat');
    warning('off', 'Octave:variable-switch-label');
    warning('off', 'Octave:fortran-indexing');
    warning('off', 'Octave:classdef-to-struct');
    warning('off', 'Octave:legacy-function'); % For strmatch and isdir

    % The following is necessary because of matlab/pac-tools/+pac/+bgp/{get,set}.m
    % which triggers this bug: https://savannah.gnu.org/bugs/?46849
    % The warning can be reenabled once the bug is fixed (and the files
    % mentioned above can then be moved back to matlab/+pac/+bgp/,
    % see #1817).
    warning('off', 'Octave:shadowed-function');
else
    % In MATLAB >= 7.7, don't display a warning if we use deprecated
    % interface to set seed of random number generators
    warning('off', 'MATLAB:RandStream:ActivatingLegacyGenerators');
end

lastwarn('') % Reset lastwarn, used to determine later if a warning has been emitted
