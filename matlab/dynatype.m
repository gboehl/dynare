function dynatype (s,var_list)
% function dynatype (s,var_list)
% This command saves the simulation results in a text file. The name of each
% variable preceeds the corresponding results. This command must follow SIMUL.
%
% INPUTS
%   s:         filename
%   var_list:  vector of selected endogenous variables
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2020 Dynare Team
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

global M_ oo_

fid = fopen(s, 'w');

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

for i = 1:length(var_list)
    idx = strcmp(var_list{i}, M_.endo_names);
    if any(idx)
        fprintf(fid, '%s\n', M_.endo_names{idx});
        fprintf(fid, '%15.8g\n', oo_.endo_simul(idx,:)');
    else
        idx = strcmp(var_list{i}, M_.exo_names);
        if any(idx)
            fprintf(fid, '%s\n', M_.exo_names{idx});
            fprintf(fid, '%15.8g\n', oo_.exo_simul(:,idx));
        else
            error(['Should not arrive here: ' var_list{i} ' not found in M_.endo_names or M_.exo_names']) ;
        end
    end
end

fclose(fid);
end
