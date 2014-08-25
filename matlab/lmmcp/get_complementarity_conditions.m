function [lb,ub,eq_index] = get_complementarity_conditions(M)

% Copyright (C) 2014 Dynare Team
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


etags = M.equations_tags;
ub = inf(M.endo_nbr,1);
lb = -ub;
eq_index = (1:M.endo_nbr)';
for i=1:size(etags,1)
    if strcmp(etags{i,2},'mcp')
        str = etags{i,3};
        kop = strfind(etags{i,3},'<');
        if ~isempty(kop)
                k = find(strcmp(strtrim(str(1:kop-1)),cellstr(M.endo_names)));
                if isempty(k)
                    error(sprintf(['Complementarity condition %s: variable %s is ' ...
                                   'not recognized',etags{i,3},b{1}]))
                end
                ub(k) = str2num(str(kop+1:end));
                eq_index(etags{i,1}) = k;
                eq_index(k) = etags{i,1};
        else
            kop = strfind(etags{i,3},'>');
            if ~isempty(kop)
                k = find(strcmp(strtrim(str(1:kop-1)),cellstr(M.endo_names)));
                if isempty(k)
                    error(sprintf(['Complementarity condition %s: variable %s is ' ...
                                   'not recognized',etags{i},b{1}]))
                end
                lb(k) = str2num(str(kop+1:end));
                eq_index(etags{i,1}) = k;
                eq_index(k) = etags{i,1};
            else
                error(sprintf(['Complementarity condition %s can''t be ' ...
                               'parsed'],etags{i,3}))
            end
        end
    end
end

