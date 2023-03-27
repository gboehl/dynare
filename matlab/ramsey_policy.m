function info = ramsey_policy(var_list)

% Copyright © 2007-2023 Dynare Team
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

global options_ oo_ M_

options_.ramsey_policy = 1;
oldoptions = options_;

[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list);

oo_.steady_state = oo_.dr.ys;

if ~options_.noprint
    disp_steady_state(M_,oo_,options_)
    for i=M_.orig_endo_nbr:M_.endo_nbr
        if strmatch('mult_', M_.endo_names{i})
            disp(sprintf('%s \t\t %g', M_.endo_names{i}, oo_.dr.ys(i)));
        end
    end
end


oo_.planner_objective_value = evaluate_planner_objective(M_,options_,oo_);

options_ = oldoptions;
