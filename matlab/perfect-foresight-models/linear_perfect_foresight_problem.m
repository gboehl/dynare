function [residuals,JJacobian] = linear_perfect_foresight_problem(y, dynamicjacobian, Y0, YT, ...
                                                  exo_simul, params, steady_state, maximum_lag, T, ny, i_cols, ...
                                                  i_cols_J1, i_cols_1, i_cols_T, i_cols_j, i_cols_0, i_cols_J0, jendo, jexog)

% Computes the residuals and the Jacobian matrix for a linear perfect foresight problem over T periods.
%
% INPUTS
% ...
%
% OUTPUTS
% ...
%
% ALGORITHM
% ...
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2015-2020 Dynare Team
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


YY = [Y0; y; YT];

residuals = zeros(T*ny,1);

z = zeros(columns(dynamicjacobian), 1);

if nargout == 2
    iJacobian = cell(T,1);
end

i_rows = 1:ny;
i_cols_J = i_cols;
offset = 0;

for it = maximum_lag+(1:T)
    z(jendo) = YY(i_cols);
    z(jexog) = transpose(exo_simul(it,:));
    residuals(i_rows) = dynamicjacobian*z;
    if nargout == 2
        if T==1 && it==maximum_lag+1
            [rows, cols, vals] = find(dynamicjacobian(:,i_cols_0));
            if size(dynamicjacobian, 1) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{1} = [rows, i_cols_J0(cols), vals];
        elseif it == maximum_lag+1
            [rows,cols,vals] = find(dynamicjacobian(:,i_cols_1));
            if size(dynamicjacobian, 1) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{1} = [offset+rows, i_cols_J1(cols), vals];
        elseif it == maximum_lag+T
            [rows,cols,vals] = find(dynamicjacobian(:,i_cols_T));
            if size(dynamicjacobian, 1) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{T} = [offset+rows, i_cols_J(i_cols_T(cols)), vals];
        else
            [rows,cols,vals] = find(dynamicjacobian(:,i_cols_j));
            if size(dynamicjacobian, 1) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{it-maximum_lag} = [offset+rows, i_cols_J(cols), vals];
            i_cols_J = i_cols_J + ny;
        end
        offset = offset + ny;
    end
    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
end

if nargout == 2
    iJacobian = cat(1,iJacobian{:});
    JJacobian = sparse(iJacobian(:,1), iJacobian(:,2), iJacobian(:,3), T*ny, T*ny);
end
