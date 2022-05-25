function message = get_error_message(info, DynareOptions)
% Returns error messages
%
% INPUTS
%   info              [double]     vector returned by resol.m
%   DynareOptions     [structure]  --> options_
% OUTPUTS
%    message          [string]     corresponding error message
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2020 Dynare Team
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

switch info(1)
    case 0
        message = '';
    case 1
        message = 'The model doesn''t determine the current variable uniquely.';
    case 2
        message = sprintf('The generalized Schur (QZ) decomposition failed. For more information, see the documentation for Lapack function dgges: info=%d, n=%d. You can also run model_diagnostics to get more information on what may cause this problem.', info(2), info(3));
    case 3
        message = 'Blanchard & Kahn conditions are not satisfied: no stable equilibrium.';
    case 4
        message = 'Blanchard & Kahn conditions are not satisfied: indeterminacy.';
    case 5
        message = 'Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.';
    case 6
        message = 'The Jacobian matrix evaluated at the steady state contains elements that are not real or are infinite.';
    case 7
        message = sprintf('One of the eigenvalues is close to 0/0 (the absolute value of numerator and denominator is smaller than %5.4f!\n If you believe that the model has a unique solution you can try to reduce the value of qz_zero_threshold.',DynareOptions.qz_zero_threshold);
    case 8
        if size(info,2)>=2 && info(2)~=0
            global M_;
            disp_string = M_.param_names{info(2)};
            for ii=1:length(info)-2
                disp_string = [disp_string, ', ', M_.param_names{info(2+ii)}];
            end
            message = ['The Jacobian contains NaNs because the following parameters are NaN: ' disp_string];
        else
            message = 'The Jacobian contains NaNs. For more information, use options_.debug.';
        end
    case 9
        message = 'k_order_pert was unable to compute the solution';
    case 10
        message = 'The Jacobian of the dynamic model contains Inf. For more information, use options_.debug.';
    case 11
        message = 'The Hessian of the dynamic model used for second order solutions must not contain Inf';
    case 12
        message = 'The Hessian of the dynamic model used for second order solutions must not contain NaN';
    case 19
        message = 'The steadystate file did not compute the steady state';
    case 20
        if DynareOptions.linear
            message = sprintf('Impossible to find the steady state (the sum of square residuals of the static equations is %5.4f). Either the model doesn''t have a steady state or there are an infinity of steady states Check whether your model is truly linear or whether there is a mistake in linearization.', info(2));
        else
            message = sprintf('Impossible to find the steady state (the sum of square residuals of the static equations is %5.4f). Either the model doesn''t have a steady state, there are an infinity of steady states, or the guess values are too far from the solution', info(2));
        end
    case 21
        message = sprintf('The steady state is complex (the sum of square residuals of imaginary parts of the steady state is %5.4f)', info(2));
    case 22
        message = 'The steady state has NaNs or Inf.';
    case 23
        message = 'Parameters have been updated in the steadystate routine and some have complex values.';
    case 24
        message = 'Parameters have been updated in the steadystate routine and some are NaNs or Inf.';
    case 25
        message = 'The solution to the static equations is not a steady state of the dynamic model: verify that the equations tagged by [static] and [dynamic] are consistent';
    case 26
        message = 'The loglinearization of the model cannot be performed, because the steady state is not strictly positive.';
    case 30
        message = 'Ergodic variance can''t be computed.';
    case 41
        message = 'one (many) parameter(s) do(es) not satisfy the lower bound';
    case 42
        message = 'one (many) parameter(s) do(es) not satisfy the upper bound';
    case 43
        message = 'Covariance matrix of structural shocks is not positive definite';
    case 44 %DsgeLikelihood_hh / dsge_likelihood
        message = 'The covariance matrix of the measurement errors is not positive definite.';
    case 45 %DsgeLikelihood_hh / dsge_likelihood
        message = 'Likelihood is not a number (NaN) or a complex number';
    case 46 %DsgeLikelihood_hh / dsge_likelihood
        message = 'Likelihood is a complex number';
    case 47 %DsgeLikelihood_hh / dsge_likelihood
        message = 'Prior density is not a number (NaN)';
    case 48 %DsgeLikelihood_hh / dsge_likelihood
        message = 'Prior density is a complex number';
    case 49
        message = 'The model violates one (many) endogenous prior restriction(s)';
    case 50
        message = 'Likelihood is Inf';
    case 51
        message = sprintf('\n The dsge_prior_weight is dsge_var=%5.4f, but must be at least %5.4f for the prior to be proper.\n You are estimating a DSGE-VAR model, but the value of the dsge prior weight is too low!', info(2), info(3));
    case 52 %dsge_var_likelihood
        message = 'You are estimating a DSGE-VAR model, but the implied covariance matrix of the VAR''s innovations, based on artificial and actual sample is not positive definite!';
    case 53 %dsge_var_likelihood
        message = 'You are estimating a DSGE-VAR model, but the implied covariance matrix of the VAR''s innovations, based on the artificial sample, is not positive definite!';
    case 55
        message = 'Fast Kalman filter only works with stationary models [lik_init=1] or stationary observables for non-stationary models [lik_init=3]';
    case 61 %Discretionary policy
        message = 'Discretionary policy: maximum number of iterations has been reached. Procedure failed.';
    case 62
        message = 'Discretionary policy: some eigenvalues greater than options_.qz_criterium. Model potentially unstable.';
    case 63
        message = 'Discretionary policy: NaN elements are present in the solution. Procedure failed.';
    case 64
        message = 'discretionary_policy: the derivatives of the objective function contain NaN.';
    case 65
        message = 'discretionary_policy: the model must be written in deviation form and not have constant terms or an analytical steady state meeds to be provided.';
    case 66
        message = 'discretionary_policy: the objective function must have zero first order derivatives.';
    case 71
        message = 'Calibrated covariance of the structural errors implies correlation larger than +-1.';
    case 72
        message = 'Calibrated covariance of the measurement errors implies correlation larger than +-1.';
        % Aim Code Conversions by convertAimCodeToInfo.m
    case 81
        message = ['Ramsey: The solution to the static first order conditions for optimal policy could not be found. Either the model' ...
            ' doesn''t have a steady state, there are an infinity of steady states, ' ...
            ' or the guess values are too far from the solution'];
    case 82
        message = 'Ramsey: The steady state computation resulted in NaN in the static first order conditions for optimal policy';
    case 83
        message = 'Ramsey: The steady state computation resulted in NaN in the auxiliary equations for optimal policy';
    case 84
        message = 'Ramsey: The steady state file computation for the Ramsey problem resulted in NaNs at the initial values of the instruments';
    case 85
        message = 'Ramsey: The steady state file does not solve the static first order conditions conditional on the instruments.';
    case 86
        message = 'Ramsey: The steady state file provides complex numbers conditional on the instruments.';
    case 87
        message = 'Ramsey: The maximum number of iterations has been reached. Try increasing maxit.';
    case 102
        message = 'Aim: roots not correctly computed by real_schur';
    case 103
        message = 'Aim: too many explosive roots: no stable equilibrium';
    case 135
        message = 'Aim: too many explosive roots, and q(:,right) is singular';
    case 104
        message = 'Aim: too few explosive roots: indeterminacy';
    case 145
        message = 'Aim: too few explosive roots, and q(:,right) is singular';
    case 105
        message = 'Aim: q(:,right) is singular';
    case 161
        message = 'Aim: too many exact shiftrights';
    case 162
        message = 'Aim: too many numeric shiftrights';
    case 163
        message = 'Aim: A is NAN or INF.';
    case 164
        message = 'Aim: Problem in SPEIG.';
    case 180
        message = 'SMM: simulation resulted in NaN/Inf. You may need to enable pruning.';
    case 201
        message = 'Particle Filter: Initial covariance of the states is not positive definite. Try a different nonlinear_filter_initialization';
    case 202
        message = 'Particle Filter: Initial covariance of the states based on simulation resulted in NaN/Inf. Use pruning or try a different nonlinear_filter_initialization';
    case 301
        message = 'IVF: The likelihood is Inf.';
    case 302
        message = 'IVF: The likelihood is NaN.';
    case 303
        message = 'IVF: The residuals are not 0.';        
    case 304
        message = 'IVF: The solver returned with an error code.';        
    case 305
        message = 'IVF: The returned shocks are bigger than 1e8.';        
    case 310
        message = 'Occbin: Simulation terminated with periodic solution (no convergence).';        
    case 311
        message = 'Occbin: Simulation did not converge, increase maxit or check_ahead_periods.';        
    case 312
        message = 'Occbin: Constraint(s) are binding at the end of the sample.';        
    case 320
        message = 'Piecewise linear Kalman filter: There was a problem in obtaining the likelihood.';        
    otherwise
        message = 'This case shouldn''t happen. Contact the authors of Dynare';
end