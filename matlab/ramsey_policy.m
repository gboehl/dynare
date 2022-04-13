function info = ramsey_policy(var_list)

% Copyright © 2007-2021 Dynare Team
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

%test whether specification matches
inst_nbr = size(options_.instruments,1);
if inst_nbr~=0
    orig_endo_aux_nbr = M_.orig_endo_nbr + min(find([M_.aux_vars.type] == 6)) - 1;
    implied_inst_nbr = orig_endo_aux_nbr - M_.orig_eq_nbr;
    if inst_nbr>implied_inst_nbr
        error('You have specified more instruments than there are omitted equations')
    elseif inst_nbr<implied_inst_nbr
        error('You have specified fewer instruments than there are omitted equations')
    end
else
    if options_.steadystate_flag
        error('You have specified a steady state file, but not provided an instrument. Either delete the steady state file or provide an instrument')
    end
end

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
