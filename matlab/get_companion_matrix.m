function get_companion_matrix(var_model_name, pac_model_name)
%function get_companion_matrix(var_model_name)

% Gets the companion matrix associated with the var specified by
% var_model_name. Output stored in cellarray oo_.var.(var_model_name).H.
%
% INPUTS
% - var_model_name   [string]        the name of the VAR model
%
% OUTPUTS
% - None

% Copyright (C) 2018 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global oo_ M_

get_ar_ec_matrices(var_model_name);

% Get the number of lags
p = size(oo_.var.(var_model_name).ar, 3);

% FIXME
while all(oo_.var.(var_model_name).ar(:,:,p)==0)
    p = p-1;
    oo_.var.(var_model_name).ar = oo_.var.(var_model_name).ar(:,:,1:p);
end

% Get the number of variables
n = length(oo_.var.(var_model_name).ar(:,:,1));

if all(~oo_.var.(var_model_name).ec(:))
    % The auxiliary model is a VAR model.
    M_.pac.(pac_model_name).auxmodel = 'var';
    if isempty(M_.pac.(pac_model_name).undiff_eqtags)
        % Build the companion matrix (standard VAR)
        oo_.var.(var_model_name).CompanionMatrix = zeros(n*p);
        oo_.var.(var_model_name).CompanionMatrix(1:n,1:n) = oo_.var.(var_model_name).ar(:,:,1);
        if p>1
            for i=2:p
                oo_.var.(var_model_name).CompanionMatrix(1:n,(i-1)*n+(1:n)) = oo_.var.(var_model_name).ar(:,:,i);
                oo_.var.(var_model_name).CompanionMatrix((i-1)*n+(1:n),(i-2)*n+(1:n)) = eye(n);
            end
        end
    else
        error('You should not use undiff option in this model!')
    end
else
    % The auxiliary model is a VECM model.
    M_.pac.(pac_model_name).auxmodel = 'vecm';
    if ~isempty(M_.pac.(pac_model_name).undiff_eqtags)
        % REMARK It is assumed that the equations with undiff option are the
        %        ECM equations. By complementarity, the other equations are
        %        the trends appearing in the error correction terms. We
        %        assume that the model can be cast in the following form:
        %
        %        Δ X_t = A_0 (X_{t-1} - Z_{t-1}) + Σ_{i=1}^p A_i Δ X_{t-i} + ϵ_t
        %
        %        Z_t = Z_{t-1} + η_t
        %
        %        We first recast the equation into this representation, and
        %        we rewrite the model in levels (we integrate the first set
        %        of equations) to rewrite the model as a VAR(1) model. Let
        %        Y_t = [X_t; Z_t] be the vertical concatenation of vectors
        %        X_t (variables with EC) and Z_t (trends). We have
        %
        %        Y_t = Σ_{i=1}^{p+1} B_i Y_{t-i} + [ε_t; η_t]
        %
        %        with
        %
        %               B_1 = [I+Λ+A_1, -Λ; 0, I]
        %
        %               B_i = [A_i-A_{i-1}, 0; 0, 0]   for i = 2,..., p
        %        and
        %               B_{p+1} = -[A_p, 0; 0, 0]
        %
        %        where the dimensions of I and 0 matrices can easily be
        %        deduced from the number of EC and trend equations.
        %
        % Get the indices of the equations with error correction terms.
        m = length(M_.pac.(pac_model_name).undiff_eqtags);
        q = length(M_.var.(var_model_name).eqn)-m;
        ecm_eqnums = zeros(m, 1);
        for i=1:m
            number = get_equation_number_by_tag(M_.pac.(pac_model_name).undiff_eqtags{i});
            if number>0
                ecm_eqnums(i) = number;
            else
                error('%s is not declared as an equation in the model block!', M_.pac.(pac_model_name).undiff_eqtag{i})
            end
        end
        % Check that the lhs of candidate ecm equations are at least first differences.
        difference_orders_in_error_correction_eq = zeros(m, 1);
        for i=1:m
            difference_orders_in_error_correction_eq(i) = get_difference_order(M_.var.(var_model_name).lhs(ecm_eqnums(i)));
        end
        if any(~difference_orders_in_error_correction_eq)
            error('Model %s is not a VECM model! LHS variables should be in difference', var_model_name)
        end
        % Get the indices of the trend equations.
        trend_eqnums = transpose(setdiff(M_.var.(var_model_name).eqn, ecm_eqnums));
        % Get the trend variables indices (lhs variables in trend equations).
        [id, id_trend_in_var, id2] = intersect(M_.var.(var_model_name).eqn, trend_eqnums);
        trend_variables = M_.var.(var_model_name).lhs(id_trend_in_var);
        % Get the rhs variables in trend equations.
        trend_autoregressive_variables = zeros(q, 1);
        for i=1:q
            % Check that there is only one variable on the rhs and update trend_autoregressive_variables.
            v = M_.var.(var_model_name).rhs.vars_at_eq{id_trend_in_var(i)}.var;
            if ~(length(v)==1)
                error('A trend equation (%s) must have only one variable on the RHS!', M_.var.(var_model_name).eqtags{trend_eqnums(i)})
            end
            trend_autoregressive_variables(i) = v;
            % Check that the variables on lhs and rhs have the same difference orders.
            if get_difference_order(trend_variables(i))~=get_difference_order(trend_autoregressive_variables(i))
                error('In a trend equation (%s) LHS and RHS variables must have the same difference orders!', M_.var.(var_model_name).eqtags{trend_eqnums(i)})
            end
            % Check that the trend equation is autoregressive.
            if isdiff(v)
                if ~M_.aux_vars(get_aux_variable_id(v)).type==9
                    error('In a trend equation (%s) RHS variable must be lagged LHS variable!', M_.var.(var_model_name).eqtags{trend_eqnums(i)})
                else
                    if M_.aux_vars(get_aux_variable_id(v)).orig_index~=trend_variables(i)
                        error('In a trend equation (%s) RHS variable must be lagged LHS variable!', M_.var.(var_model_name).eqtags{trend_eqnums(i)})
                    end
                end
            else
                if get_aux_variable_id(v) && M_.aux_vars(get_aux_variable_id(v)).endo_index~=trend_variables(i)
                    error('In a trend equation (%s) RHS variable must be lagged LHS variable!', M_.var.(var_model_name).eqtags{trend_eqnums(i)})
                end
            end
        end
        ecm_trend_eq = trend_eqnums;
        % Get the EC matrix (the EC term is assumend to be in t-1).
        %
        % TODO: Check that the EC term is the difference between the
        %       endogenous variable and the trend variable.
        %
        A0 = oo_.var.(var_model_name).ec(ecm_eqnums,:,1);
        % Get the AR matrices.
        AR = oo_.var.(var_model_name).ar(ecm_eqnums,ecm_eqnums,:);
        % Build B matrices (VAR in levels)
        B = zeros(n, n, p+1);
        B(ecm_eqnums,ecm_eqnums,1) = eye(m)+A0+AR(:,:,1);
        B(ecm_eqnums,ecm_trend_eq) = -A0;
        B(ecm_trend_eq,ecm_trend_eq) = eye(q);
        for i=2:p
            B(ecm_eqnums,ecm_eqnums,i) = AR(:,:,i)-AR(:,:,i-1);
        end
        B(ecm_eqnums,ecm_eqnums,p+1) = -AR(:,:,p);
        % Write Companion matrix
        oo_.var.(var_model_name).CompanionMatrix = zeros(size(B, 1)*size(B, 3));
        for i=1:p
            oo_.var.(var_model_name).CompanionMatrix(1:n, (i-1)*n+(1:n)) = B(:,:,i);
            oo_.var.(var_model_name).CompanionMatrix(i*n+(1:n),(i-1)*n+(1:n)) = eye(n);
        end
        oo_.var.(var_model_name).CompanionMatrix(1:n, p*n+(1:n)) = B(:,:,p+1);
    else
        error('It is not possible to cast the VECM model in a companion representation! Use undiff option.')
    end
end