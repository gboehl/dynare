function ys1 = add_auxiliary_variables_to_steadystate(ys,aux_vars,fname, ...
                                                  exo_steady_state, exo_det_steady_state,params, byte_code)
% Add auxiliary variables to the steady state vector

% Copyright © 2009-2023 Dynare Team
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

global M_ options_

n = length(aux_vars);
ys1 = [ys;zeros(n,1)];

for i=1:n+1
    if byte_code
        res = bytecode('static','evaluate', M_, options_, ys1,...
                       [exo_steady_state; ...
                        exo_det_steady_state],params);
    else
        res = feval([fname '.static'],ys1,...
                    [exo_steady_state; ...
                     exo_det_steady_state],params);
    end
    for j=1:n
        el = aux_vars(j).endo_index;
        ys1(el) = ys1(el)-res(el);
    end
end
