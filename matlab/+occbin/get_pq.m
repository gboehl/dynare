function [p,q]=get_pq(dr)
% function [p,q]=get_pq(dr)
% Output decision rule matrices p and q in declaration order
%
% INPUTS
% - dr         [struct]     decision rules
%
% OUTPUTS
% - p           [N by N]     transition matrix ghu in declaration order
% - q           [N by N_exo] shock response matrix ghx in declaration order

% Copyright (C) 2021 Dynare Team
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

nvars = size(dr.ghx,1);
p = zeros(nvars);
p(dr.order_var,dr.state_var) = dr.ghx;
q = dr.ghu(dr.inv_order_var,:);
