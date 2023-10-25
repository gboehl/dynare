function [residuals,JJacobian] = perfect_foresight_mcp_problem(y, dynamic_function, Y0, YT, ...
                                                  exo_simul, params, steady_state, ...
                                                  maximum_lag, T, ny, i_cols, ...
                                                  i_cols_J1, i_cols_1, i_cols_T, ...
                                                  i_cols_j, i_cols_0,i_cols_J0, eq_index)
% function [residuals,JJacobian] = perfect_foresight_mcp_problem(y, dynamic_function, Y0, YT, ...
%                                            exo_simul, params, steady_state, ...
%                                            maximum_lag, T, ny, i_cols, ...
%                                            i_cols_J1, i_cols_1, i_cols_T, ...
%                                            i_cols_j,eq_index)
% Computes the residuals and the Jacobian matrix for a perfect foresight problem over T periods
% in a mixed complementarity problem context
%
% INPUTS
%   y                   [double] N*1 array, terminal conditions for the endogenous variables
%   dynamic_function    [handle] function handle to _dynamic-file
%   Y0                  [double] N*1 array, initial conditions for the endogenous variables
%   YT                  [double] N*1 array, terminal conditions for the endogenous variables
%   exo_simul           [double] nperiods*M_.exo_nbr matrix of exogenous variables (in declaration order)
%                                for all simulation periods
%   params              [double] nparams*1 array, parameter values
%   steady_state        [double] endo_nbr*1 vector of steady state values
%   maximum_lag         [scalar] maximum lag present in the model
%   T                   [scalar] number of simulation periods
%   ny                  [scalar] number of endogenous variables
%   i_cols              [double] indices of variables appearing in M_.lead_lag_incidence
%                                and that need to be passed to _dynamic-file
%   i_cols_J1           [double] indices of contemporaneous and forward looking variables
%                                appearing in M_.lead_lag_incidence
%   i_cols_1            [double] indices of contemporaneous and forward looking variables in
%                                M_.lead_lag_incidence in dynamic Jacobian (relevant in first period)
%   i_cols_T            [double] columns of dynamic Jacobian related to contemporaneous and backward-looking
%                                variables (relevant in last period)
%   i_cols_j            [double] indices of variables in M_.lead_lag_incidence
%                                in dynamic Jacobian (relevant in intermediate periods)
%   eq_index            [double] N*1 array, index vector describing residual mapping resulting
%                                from complementarity setup
% OUTPUTS
%   residuals           [double] (N*T)*1 array, residuals of the stacked problem
%   JJacobian           [double] (N*T)*(N*T) array, Jacobian of the stacked problem
% ALGORITHM
%   None
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 1996-2020 Dynare Team
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
if nargout == 2
    iJacobian = cell(T,1);
end

i_rows = 1:ny;
offset = 0;
i_cols_J = i_cols;

for it = maximum_lag+(1:T)
    if nargout == 1
        res = dynamic_function(YY(i_cols),exo_simul, params, ...
                               steady_state,it);
        residuals(i_rows) = res(eq_index);
    elseif nargout == 2
        [res,jacobian] = dynamic_function(YY(i_cols),exo_simul, params, steady_state,it);
        residuals(i_rows) = res(eq_index);
        if T==1 && it==maximum_lag+1
            [rows, cols, vals] = find(jacobian(eq_index,i_cols_0));
            if size(jacobian, 1) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{1} = [rows, i_cols_J0(cols), vals];
        elseif it == maximum_lag+1
            [rows,cols,vals] = find(jacobian(eq_index,i_cols_1));
            if numel(eq_index) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{1} = [offset+rows, i_cols_J1(cols), vals];
        elseif it == maximum_lag+T
            [rows,cols,vals] = find(jacobian(eq_index,i_cols_T));
            if numel(eq_index) == 1 % find() will return row vectors in this case
                rows = rows';
                cols = cols';
                vals = vals';
            end
            iJacobian{T} = [offset+rows, i_cols_J(i_cols_T(cols)), vals];
        else
            [rows,cols,vals] = find(jacobian(eq_index,i_cols_j));
            if numel(eq_index) == 1 % find() will return row vectors in this case
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
    JJacobian = sparse(iJacobian(:,1),iJacobian(:,2),iJacobian(:,3),T*ny,T*ny);
end
