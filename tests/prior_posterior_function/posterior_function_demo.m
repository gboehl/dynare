function output_cell =posterior_function_demo(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info)
% output_cell =posterior_function_demo(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info);
% This is an example file computing statistics on the prior/posterior draws. The
% function allows read-only access to all Dynare structures. However, those
% structures are local to this function.  Changing them will not affect
% other Dynare functions and you cannot use them to pass results to other
% Dynare functions.
% The function takes one and only one output argument: an 1 by n cell.
% Using functions like cell2mat, the contents of the cell can be easily
% transformed back to matrices. See the fs2000_posterior_function.mod for
% an example

% INPUTS
%   xparam1                      Current parameter draw
%   M_           [structure]     Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%   options_     [structure]     Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_          [structure]     Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%   estim_params_[structure]     Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   bayestopt_   [structure]     Matlab's structure describing the parameter options (initialized by dynare, see @ref{bayestopt_}).
%   dataset_     [structure]     Matlab's structure storing the dataset
%   dataset_info [structure]     Matlab's structure storing the information about the dataset

% Output
%   output_cell  [1 by n cell]   1 by n Matlab cell allowing to store any
%                                desired computation or result (strings, matrices, structures, etc.)

% Copyright © 2015-2023 Dynare Team
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


%% store the mean of the parameter draw
output_cell{1,1}=mean(xparam1);

%% compute the steady state for the parameter draw and store it
% set the parameters draws to the model structure
M_ = set_all_parameters(xparam1,estim_params_,M_);
% compute the steady state for the parameter draw written to M_
[ys,params,info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);

%set second part of output cell
output_cell{1,2}=ys';
end
