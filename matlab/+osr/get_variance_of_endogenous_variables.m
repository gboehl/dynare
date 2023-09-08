function vx1 = get_variance_of_endogenous_variables(M_,options_,dr,i_var)
% vx1 = get_variance_of_endogenous_variables(dr,i_var)
% Gets the variance of a variables subset
%
% INPUTS
%   M_                        [structure]       Dynare's model structure
%   oo_                       [structure]       Dynare's results structure
%   options_                  [structure]       Dynare's options structure
%   dr:                       [structure]       structure of decisions rules for stochastic simulations
%   i_var:                    [integer]         indices of a variables list
%
% OUTPUTS
%    vx1:                     [double]          variance-covariance matrix
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2023 Dynare Team
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

ghx = dr.ghx(i_var,:);
ghu = dr.ghu(i_var,:);
nc = size(ghx,2);
n = length(i_var);

[A,B] = kalman_transition_matrix(dr,M_.nstatic+(1:M_.nspred),1:nc);

[vx,u] = lyapunov_symm(A,B*M_.Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold, [], options_.debug);

if size(u,2) > 0
    i_stat = find(any(abs(ghx*u) < options_.schur_vec_tol,2)); %only set those variances of objective function for which variance is finite
    ghx = ghx(i_stat,:);
    ghu = ghu(i_stat,:);
else
    i_stat = (1:n)';
end

vx1 = Inf*ones(n,n);
vx1(i_stat,i_stat) = ghx*vx*ghx'+ghu*M_.Sigma_e*ghu';