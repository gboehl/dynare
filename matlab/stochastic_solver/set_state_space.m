function dr=set_state_space(dr,M_)
% dr=set_state_space(dr,M_)
% This function computes the DR ordering and inverse ordering.
% It used to compute stuff related to the state-space representation of the reduced form, hence
%  its name.

%@info:
%! @deftypefn {Function File} {[@var{dr} =} set_state_space (@var{dr},@var{M_})
%! @anchor{set_state_space}
%! @sp 1
%! Write the state space representation of the reduced form solution.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing decision and transition rules.
%! @item M_
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_})
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing decision and transition rules.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{check}, @ref{discretionary_policy_1}, @ref{dynare_estimation_init}, @ref{dyn_risky_steady_state_solver}, @ref{osr1}, @ref{partial_information/dr1_PI}, @ref{pea/pea_initialization}, @ref{stochastic_solvers}, @ref{stoch_simul}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright © 1996-2024 Dynare Team
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

max_lead = M_.maximum_endo_lead;
max_lag = M_.maximum_endo_lag;
endo_nbr = M_.endo_nbr;
lead_lag_incidence = M_.lead_lag_incidence;
klen = max_lag + max_lead + 1;

fwrd_var = find(lead_lag_incidence(max_lag+2:end,:))';
if max_lag > 0
    pred_var = find(lead_lag_incidence(1,:))';
    both_var = intersect(pred_var,fwrd_var);
    pred_var = setdiff(pred_var,both_var);
    fwrd_var = setdiff(fwrd_var,both_var);
    stat_var = setdiff([1:endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
    pred_var = [];
    both_var = [];
    stat_var = setdiff([1:endo_nbr]',fwrd_var);
end

dr.order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
dr.inv_order_var(dr.order_var) = 1:endo_nbr;

dr.transition_auxiliary_variables = [];
