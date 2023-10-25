function [endogenousvariables, exogenousvariables] = backward_model_inversion(constraints, exogenousvariables, initialconditions, endo_names, exo_names, freeinnovations, M_, options_, oo_)
% [endogenousvariables, exogenousvariables] = backward_model_inversion(constraints, exogenousvariables, initialconditions, endo_names, exo_names, freeinnovations, M_, options_, oo_)
% INPUTS
% - constraints         [dseries]        with N constrained endogenous variables from t1 to t2.
% - exogenousvariables  [dseries]        with Q exogenous variables.
% - initialconditions   [dseries]        with M endogenous variables starting before t1 (M initialcond must contain at least the state variables).
% - endo_names          [cell]           list of endogenous variable names.
% - exo_names           [cell]           list of exogenous variable names.
% - freeinstruments     [cell]           list of exogenous variable names used to control the constrained endogenous variables.
% - M_                  [structure]      describing the model
% - options_            [structure]      describing the options
% - oo_                 [structure]      storing the results
%
% OUTPUTS
% - endogenous          [dseries]
% - exogenous           [dseries]
%
% REMARKS

% Copyright © 2017-2023 Dynare Team
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

% Get indices for the free innovations.
freeinnovations_id = zeros(length(freeinnovations), 1);
if length(freeinnovations)<M_.exo_nbr
    for i=1:length(freeinnovations)
        freeinnovations_id(i) = find(strcmp(freeinnovations{i}, exo_names));
    end
else
    freeinnovations_id = transpose(1:length(exo_names));
end

nxfree = length(freeinnovations_id);

% Get indices for the the controlled and free endogenous variables.
controlledendogenousvariables_id = zeros(length(freeinnovations), 1);
if length(freeinnovations)<M_.endo_nbr
    for i=1:length(freeinnovations)
        controlledendogenousvariables_id(i) = find(strcmp(constraints.name{i}, endo_names));
    end
    freeendogenousvariables_id = setdiff(transpose(1:length(endo_names)), controlledendogenousvariables_id);
else
    controlledendogenousvariables_id = transpose(1:length(endo_names));
    freeendogenousvariables_id = [];
end

nyfree = length(freeendogenousvariables_id);
nyctrl = length(controlledendogenousvariables_id);

% Build structure to be passed to the objective function.
ModelInversion.nyfree = nyfree;
ModelInversion.nyctrl = nyctrl;
ModelInversion.nxfree = nxfree;
ModelInversion.y_constrained_id = controlledendogenousvariables_id;
ModelInversion.y_free_id = freeendogenousvariables_id;
ModelInversion.x_free_id = freeinnovations_id;
ModelInversion.J_id = [M_.endo_nbr+ModelInversion.y_free_id ; 3*M_.endo_nbr+ModelInversion.x_free_id];

% Get function handles to the dynamic model routines.
dynamic_resid = str2func([M_.fname '.sparse.dynamic_resid']);
dynamic_g1 = str2func([M_.fname '.sparse.dynamic_g1']);

% Initialization of the returned simulations (endogenous variables).
Y = NaN(M_.endo_nbr, nobs(constraints));
Y = [transpose(initialconditions{M_.endo_names{:}}(constraints.dates(1)-1).data), Y];
for i=1:nyctrl
    Y(controlledendogenousvariables_id(i),2:end) = transpose(constraints.data(:,i));
end

% Exogenous variables.
X = exogenousvariables{exo_names{:}}.data;

% Inversion of the model, solvers for the free endogenous and exogenous variables (call a Newton-like algorithm in each period).
ity = 2;
itx = find(exogenousvariables.dates==constraints.dates(1));
for t = 1:nobs(constraints)
    % Set the lagged values of the endogenous variables.
    ylag = Y(:,ity-1);
    % Set the current values of the constrained endogenous variables.
    ycur = Y(controlledendogenousvariables_id,ity);
    % Vector z gather the free endogenous variables (initialized with lagged
    % values) and the free exogenous variables (initialized with 0).
    z = [Y(freeendogenousvariables_id,ity-1); zeros(nxfree, 1)];
    % Solves for z.
    [z, failed, ~, ~, errorcode] = dynare_solve(@dynamic_backward_model_for_inversion, z, ...
                                                options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                                options_, dynamic_resid, dynamic_g1, ylag, ycur, X(itx,:), M_.params, oo_.steady_state, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr, ModelInversion);
    if failed
        error('Nonlinear solver failed with errorcode=%i', errorcode)
    end
    % Update the matrix of exogenous variables.
    X(itx,freeinnovations_id) = z(nyfree+(1:nxfree));
    % Update the matrix of endogenous variables.
    if nyfree
        Y(freeendogenousvariables_id,ity) = z(1:nyfree);
    end
    % Increment counters
    ity = ity+1;
    itx = itx+1;
end

endogenousvariables = dseries(Y(:,2:end)', constraints.dates(1), endo_names);
exogenousvariables = dseries(X(find(exogenousvariables.dates==constraints.dates(1))+(0:(nobs(constraints)-1)),:), constraints.dates(1), exo_names);


function [r, J] = dynamic_backward_model_for_inversion(z, dynamic_resid, dynamic_g1, ylag, ycur, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, ModelInversion)

endo_nbr = ModelInversion.nyfree+ModelInversion.nyctrl;

% Set up y
y = [ylag; zeros(2*endo_nbr, 1)];

y(endo_nbr+ModelInversion.y_constrained_id) = ycur;
if ModelInversion.nyfree
    y(endo_nbr+ModelInversion.y_free_id) = z(1:ModelInversion.nyfree);
end

% Update x
x(ModelInversion.x_free_id) = z(ModelInversion.nyfree+(1:ModelInversion.nxfree));

[r, T_order, T] = dynamic_resid(y, x, params, steady_state);

if nargout>1
    Jacobian = dynamic_g1(y, x, params, steady_state, sparse_rowval, ...
                          sparse_colval, sparse_colptr, T_order, T);
    J = Jacobian(:, ModelInversion.J_id);
end
