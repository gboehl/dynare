function [r, J] = dynamic_forward_model_for_simulation(z, dynamicmodel, ylead, x, params, steady_state, it_)

% Dynamic routine's wrapper used by dynare_solve.

% Copyright © 2022 Dynare Team
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

if nargout>1
        % Compute residuals and jacobian of the full dynamic model.
    [r, J] = feval(dynamicmodel, [z; ylead], x, params, steady_state, it_);
    J = J(:,1:length(z)); % Remove derivatives with respect to shocks.
else
    % Compute residuals.
    r = feval(dynamicmodel, [z; ylead], x, params, steady_state, it_);
end