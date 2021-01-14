function [out1, out2] = method_of_moments_objective_function_gradient_helper(xparam1, Bounds, oo_, estim_params_, M_, options_mom_)
% [out1, out2] = method_of_moments_objective_function_gradient_helper(xparam1, Bounds, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This helper function evaluates the objective function for GMM/SMM estimation and
% outputs the function value fval at first and the gradient df at second place,
% needed for gradient-based optimizers if analytic_jacobian option is set
% =========================================================================
% INPUTS
%   o xparam1:                  current value of estimated parameters as returned by set_prior()
%   o Bounds:                   structure containing parameter bounds
%   o oo_:                      structure for results
%   o estim_params_:            structure describing the estimated_parameters
%   o M_                        structure describing the model
%   o options_mom_:             structure information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS: dependent on the optimizer calling this function, see below
% -------------------------------------------------------------------------
% This function calls
%  o method_of_moments_objective_function.m
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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
% -------------------------------------------------------------------------
% This function is called by
%  o method_of_moments.m
%  o dynare_minimize_objective.m
% -------------------------------------------------------------------------
% This function calls
%  o method_of_moments_objective_function
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% =========================================================================
[fval, info, exit_flag, junk1, junk2, oo_, M_, options_mom_, df] = method_of_moments_objective_function(xparam1, Bounds, oo_, estim_params_, M_, options_mom_);

switch options_mom_.current_optimizer
    case 1 %fmincon
        out1=fval; out2=df;
    case 3 %fminunc
        out1=fval; out2=df;
    case 13 %lsqnonlin
        out1=fval; out2=df;
    otherwise
        error('Method of Moments: analytic_jacobian option is currently only supported by mode_compute 1, 3 and 13');
end

end%main function end

