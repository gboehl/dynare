function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M,options)

% function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M,options)
% Evaluates the static model
%
% INPUTS
%   ys                        vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state
%   params                    vector           parameters
%   M                         struct           model structure
%   options                   struct           options
%
% OUTPUTS
%   residuals                 vector           residuals when ys is not
%                                              the steady state
%   check1                    scalar           error flag
%   jacob                     matrix           Jacobian of static model
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2021 Dynare Team
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
if options.bytecode
    if nargout<3
        [residuals]= bytecode('evaluate','static',ys,...
                         exo_ss, params, ys, 1);
    else
        [residuals, junk]= bytecode('evaluate','static',ys,...
            exo_ss, params, ys, 1);
        jacob = junk.g1;
    end      
else
    fh_static = str2func([M.fname '.static']);
    if options.block
        residuals = zeros(M.endo_nbr,1);
        T = NaN(M.block_structure_stat.tmp_nbr, 1);
        for b = 1:length(M.block_structure_stat.block)
            [r, yy, T] = feval(fh_static,b,ys,exo_ss,params,T);
            if M.block_structure_stat.block(b).Simulation_Type == 1 || ... % evaluateForward
               M.block_structure_stat.block(b).Simulation_Type == 2        % evaluateBackward
                vidx = M.block_structure_stat.block(b).variable;
                r = yy(vidx) - ys(vidx);
            end
            residuals(M.block_structure_stat.block(b).equation) = r;
        end
        if nargout==3
            jacob=NaN(length(ys));
        end
    else
        if nargout<3
            residuals = feval(fh_static,ys,exo_ss,params);
        else
            [residuals, jacob] = feval(fh_static,ys,exo_ss,params);
        end
    end
end
