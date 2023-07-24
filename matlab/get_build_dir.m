function p = get_build_dir(dynareroot)
% Returns a Meson build directory if one is found.
% Otherwise returns an empty value.

% Copyright © 2023 Dynare Team
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

envvar = getenv('DYNARE_BUILD_DIR');
default_matlab = [ dynareroot '..' filesep 'build-matlab' ];
default_octave = [ dynareroot '..' filesep 'build-octave' ];

if ~isempty(envvar)
    p = envvar;
elseif ~isoctave && exist(default_matlab, 'dir')
    p = default_matlab;
elseif isoctave && exist(default_octave, 'dir')
    p = default_octave;
else
    p = [];
end

end
