function [loss,info,exit_flag,vx,junk]=objective(x,M_, oo_, options_,i_params,i_var,weights)
% objective function for optimal simple rules (OSR)
% INPUTS
%   x                         vector           values of the parameters
%                                              over which to optimize
%   M_                        [structure]      Dynare's model structure
%   oo_                       [structure]      Dynare's results structure
%   options_                  [structure]      Dynare's options structure
%   i_params                  vector           index of optimizing parameters in M_.params
%   i_var                     vector           variables indices
%   weights                   vector           weights in the OSRs
%
% OUTPUTS
%   loss                      scalar           loss function returned to solver
%   info                      vector           info vector returned by resol
%   exit_flag                 scalar           exit flag returned to solver
%   vx                        vector           variances of the endogenous variables
%   junk                      empty            dummy output for conformable
%                                              header
%
% SPECIAL REQUIREMENTS
%   none
% Copyright © 2005-2023 Dynare Team
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

junk = [];
exit_flag = 1;
vx = [];
% set parameters of the policiy rule
M_.params(i_params) = x;

% don't change below until the part where the loss function is computed
[dr,info] = resol(0,M_,options_,oo_);

if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85
        loss = 1e8;
        info(4)=info(2);
        return
    else
        loss = 1e8;
        info(4)=0.1;
        return
    end
end

vx = osr.get_variance_of_endogenous_variables(M_,options_,dr,i_var);
loss = full(weights(:)'*vx(:));