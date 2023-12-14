function clear_persistent_variables(folder, writelistofroutinestobecleared)
% clear_persistent_variables(folder, writelistofroutinestobecleared)
% Clear all the functions with persistent variables in directory folder (and subdirectories).

% Copyright © 2015-2019 Dynare Team
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
    writelistofroutinestobecleared = false;
end

if nargin<1 || isempty(folder)
    folder = pwd();
end

DYNARE_FOLDER = strrep(which('dynare'),'dynare.m','');

if writelistofroutinestobecleared
    if ~exist('list_of_functions_to_be_cleared.m') || isolder(sprintf('%slist_of_functions_to_be_cleared.m', DYNARE_FOLDER), DYNARE_FOLDER)
        if isunix() || ismac()
            [~, output] = system(sprintf('grep -lr ^persistent %s', folder));
            list_of_files = strsplit(output);
            list_of_files(find(cellfun(@isempty, list_of_files))) = [];
        else
            [~, output] = system(sprintf('findstr /B/S/M persistent %s\\*', folder));
            list_of_files = strsplit(output);
            list_of_files(find(cellfun(@isempty, list_of_files))) = [];
            i = 1; mobius = true;
            while mobius
                if i>length(list_of_files)
                    break
                end
                if ismember(list_of_files{i},{'FINDSTR:', 'ignored', '//'})
                    list_of_files(i) = [];
                else
                    i = i + 1;
                end
            end
        end
        [~, list_of_functions, ~] = cellfun(@fileparts, list_of_files, 'UniformOutput',false);
        cellofchar2mfile(sprintf('%slist_of_functions_to_be_cleared.m', DYNARE_FOLDER), list_of_functions)
    end
    return
end

list_of_functions_to_be_cleared;
clear(list_of_functions{:});

function cellofchar2mfile(fname, c, cname)
% Write a cell of char in a matlab script.
%
% INPUTS
% - fname [string] name of the file where c is to be saved.
% - c     [cell]   a two dimensional cell of char.
%
% OUTPUTS
% None.



[pathstr,name,ext] = fileparts(fname);

if isempty(ext)
    fname = [pathstr, name, '.m'];
else
    if ~isequal(ext, '.m')
        error(['The first argument needs to be the name of a matlab script (with an .m extension)!'])
    end
end

if ~iscell(c)
    error('The second input argument must be a cell!')
end

if ndims(c)>2
    error(['The cell passed has a second argument cannot have more than two dimensions!'])
end

variablename = inputname(2);

if isempty(variablename) && nargin<3
    error(['You must pass the name of the cell (second input argument) as a string in the third input argument!'])
end

if nargin>2
    if isvarname(cname)
        variablename = cname;
    else
        error('The third input argument must be a valid variable name!')
    end
end

fid = fopen(fname,'w');
fprintf(fid, '%s = %s;', variablename, writecellofchar(c));
fclose(fid);
