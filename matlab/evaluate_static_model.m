function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M_,options_)

% function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M_,options_)
% Evaluates the static model
%
% INPUTS
%   ys                        vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state
%   params                    vector           parameters
%   M_                        struct           model structure
%   options_                  struct           options
%
% OUTPUTS
%   residuals                 vector           residuals when ys is not
%                                              the steady state
%   check1                    scalar           error flag
%   jacob                     matrix           Jacobian of static model
%
% SPECIAL REQUIREMENTS
%   none

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

check1 = 0;
if options_.bytecode
    if nargout<3
        [residuals]= bytecode('evaluate', 'static', M_, options_, ys, ...
                         exo_ss, params, ys, 1);
    else
        [residuals, junk]= bytecode('evaluate', 'static', M_, options_, ys, ...
            exo_ss, params, ys, 1);
        jacob = junk.g1;
    end      
else
    [residuals, T_order, T] = feval([M_.fname '.sparse.static_resid'], ys, exo_ss, params);
    if nargout >= 3
        jacob = feval([M_.fname '.sparse.static_g1'], ys, exo_ss, params, M_.static_g1_sparse_rowval, M_.static_g1_sparse_colval, M_.static_g1_sparse_colptr, T_order, T);
    end
end
