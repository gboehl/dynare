function var_decomp=compute_variance_decomposition(M_,options_,var_stationary,A,A_stationary,ghu,ghu_states_only,stationary_vars,index_stationary_vars,nvars_requested)
% function var_decomp=compute_variance_decomposition(M_,options_,var_stationary,A,A_stationary,ghu,ghu_states_only,stationary_vars,index_stationary_vars,nvars_requested)
% Computes the theoretical auto-covariances, Gamma_y, for an AR(p) process
% with coefficients dr.ghx and dr.ghu and shock variances Sigma_e
% for a subset of variables ivar.
% Theoretical HP-filtering and band-pass filtering is available as an option
%
% INPUTS
%   M_                      [structure]     Global dynare's structure, description of the DSGE model.
%   options_                [structure]     Global dynare's structure.
%   var_stationary          [double]        unconditional variance of stationary
%                                           variables
%   A                       [double]        State transition matrix
%   A_stationary            [double]        Reaction of stationary variables to states
%   ghu                     [double]        Reaction of variables to shocks
%   ghu_states_only         [double]        Reaction of stationary variables to shocks
%   stationary_vars         [integer]       index of stationary vars in requested output
%   index_stationary_vars   [integer]       index of stationary vars in decision rules
%   nvars_requested         [integer]       number of originally requested variables
%
% OUTPUTS
%   stationary_vars   [double]      [#stationary vars by shocks] Matrix containing the variance decomposition 
%
% Copyright © 2001-2023 Dynare Team
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

if M_.exo_nbr == 1
    var_decomp = ones(size(ghu,1),1);
else
    var_decomp = NaN(nvars_requested,M_.exo_nbr);
    cs = get_lower_cholesky_covariance(M_.Sigma_e,options_.add_tiny_number_to_cholesky);
    b1 = ghu_states_only*cs;
    b2 = ghu(index_stationary_vars,:)*cs;
    variance_sum_loop = 0;
    for i=1:M_.exo_nbr
        if i==1 % first time, do Schur decomposition
            variance_states = lyapunov_symm(A,b1(:,i)*b1(:,i)',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,1,options_.debug);
        else
            variance_states = lyapunov_symm(A,b1(:,i)*b1(:,i)',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,2,options_.debug);
        end        
        vx2 = diag(A_stationary*variance_states*A_stationary'+b2(:,i)*b2(:,i)');
        var_decomp(stationary_vars,i) = vx2;
        variance_sum_loop = variance_sum_loop +vx2; %track overall variance over shocks
    end
    if ~options_.pruning && max(abs(variance_sum_loop-var_stationary)./var_stationary) > 1e-4 && max(abs(variance_sum_loop-var_stationary))>1e-7
        warning(['Aggregate variance and sum of variances by shocks ' ...
            'differ by more than 0.01 %'])
    end
    var_decomp(stationary_vars,:) = var_decomp(stationary_vars,:)./repmat(variance_sum_loop,1,M_.exo_nbr);
end