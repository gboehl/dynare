function [A0, AR, B] = get_companion_matrix_legacy(auxiliary_model_name, auxiliary_model_type)

% Gets the companion VAR representation of a PAC auxiliary model.
% Depending on the nature of this auxiliary model the output is
% saved in oo_.{var,trend_component}.(auxiliary_model_name).CompanionMatrix
%
% INPUTS
% - auxiliary_model_name     [string]    the name of the auxiliary model
% - auxiliary_model_type     [string]    the type of the auxiliary model
%                                        ('var' or 'trend_component')
%
% OUTPUTS
% - None

% Copyright © 2018 Dynare Team
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

global oo_ M_

if nargin<2
    if isfield(M_, 'var') && isfield(M_.var, auxiliary_model_name)
        auxiliary_model_type = 'var';
    elseif isfield(M_, 'trend_component') && isfield(M_.trend_component, auxiliary_model_name)
        auxiliary_model_type = 'trend_component';
    else
        error('Unknown type of auxiliary model.')
    end
end

if nargout
    A0 = [];
    AR = [];
    B = [];
end

get_ar_ec_matrices(auxiliary_model_name, auxiliary_model_type);

% Get the number of lags
p = size(oo_.(auxiliary_model_type).(auxiliary_model_name).ar, 3);

% Get the number of variables
n = length(oo_.(auxiliary_model_type).(auxiliary_model_name).ar(:,:,1));

switch auxiliary_model_type
  case 'var'
    oo_.var.(auxiliary_model_name).CompanionMatrix = zeros(n*p);
    oo_.var.(auxiliary_model_name).CompanionMatrix(1:n,1:n) = oo_.var.(auxiliary_model_name).ar(:,:,1);
    for i=2:p
        oo_.var.(auxiliary_model_name).CompanionMatrix(1:n,(i-1)*n+(1:n)) = oo_.var.(auxiliary_model_name).ar(:,:,i);
        oo_.var.(auxiliary_model_name).CompanionMatrix((i-1)*n+(1:n),(i-2)*n+(1:n)) = eye(n);
    end
    AR = oo_.var.(auxiliary_model_name).ar;
    M_.var.(auxiliary_model_name).list_of_variables_in_companion_var = M_.endo_names(M_.var.(auxiliary_model_name).lhs);
  case 'trend_component'
    % Get number of trends.
    q = sum(M_.trend_component.(auxiliary_model_name).targets);
    % Get the number of equations with error correction.
    m  = n-q;
    % Get the indices of trend and EC equations in the auxiliary model.
    target_eqnums_in_auxiliary_model = find(M_.trend_component.(auxiliary_model_name).targets);
    ecm_eqnums_in_auxiliary_model = find(~M_.trend_component.(auxiliary_model_name).targets);
    % REMARK It is assumed that the non trend equations are the error correction
    %        equations. We assume that the model can be cast in the following form:
    %
    %        Δ Xₜ₋₁ = A₀ (Xₜ₋₁ - Zₜ₋₁) + Σᵢ₌₁ᵖ Aᵢ Δ Xₜ₋ᵢ + ϵₜ
    %
    %        Zₜ = Zₜ₋₁ + ηₜ
    %
    %        We first recast the equation into this representation, and
    %        we rewrite the model in levels (we integrate the first set
    %        of equations) to rewrite the model as a VAR(1) model. Let
    %        Yₜ = [Xₜ; Zₜ] be the vertical concatenation of vectors
    %        Xₜ (variables with EC) and Zₜ (trends). We have
    %
    %        Yₜ = Σᵢ₌₁ᵖ⁺¹ Bᵢ Yₜ₋ᵢ + [εₜ; ηₜ]
    %
    %        with
    %
    %               B₁ = [I+Λ+A₁, -Λ; 0, I]
    %
    %               Bᵢ = [Aᵢ-Aᵢ₋₁, 0; 0, 0]   for i = 2,…, p
    %        and
    %               Bₚ₊₁ = -[Aₚ, 0; 0, 0]
    %
    %        where the dimensions of I and 0 matrices can easily be
    %        deduced from the number of EC and trend equations.
    % Check that the lhs of candidate ecm equations are at least first differences.
    for i=1:m
        if ~get_difference_order(M_.trend_component.(auxiliary_model_name).lhs(ecm_eqnums_in_auxiliary_model(i)))
            error('Model %s is not a Trend component  model! LHS variables should be in difference', auxiliary_model_name)
        end
    end
    % Reorder target_eqnums_in_auxiliary_model to ensure that the order of
    % the trend variables matches the order of the error correction
    % variables.
    [~,reorder] = ismember(M_.trend_component.(auxiliary_model_name).lhs(target_eqnums_in_auxiliary_model), ...
                           M_.trend_component.(auxiliary_model_name).target_vars(find(M_.trend_component.(auxiliary_model_name).target_vars>0)));
    target_eqnums_in_auxiliary_model = target_eqnums_in_auxiliary_model(reorder);
    % Get the EC matrix (the EC term is assumend to be in t-1).
    %
    % TODO: Check that the EC term is the difference between the
    %       endogenous variable and the trend variable.
    %
    A0 = oo_.trend_component.(auxiliary_model_name).ec(ecm_eqnums_in_auxiliary_model,:,1);
    % Get the AR matrices.
    AR = oo_.trend_component.(auxiliary_model_name).ar(ecm_eqnums_in_auxiliary_model,ecm_eqnums_in_auxiliary_model,:);
    % Build B matrices (VAR in levels)
    B(ecm_eqnums_in_auxiliary_model,ecm_eqnums_in_auxiliary_model,1) = eye(m)+A0+AR(:,:,1);
    B(ecm_eqnums_in_auxiliary_model,target_eqnums_in_auxiliary_model) = -A0;
    B(target_eqnums_in_auxiliary_model,target_eqnums_in_auxiliary_model) = eye(q);
    for i=2:p
        B(ecm_eqnums_in_auxiliary_model,ecm_eqnums_in_auxiliary_model,i) = AR(:,:,i)-AR(:,:,i-1);
    end
    B(ecm_eqnums_in_auxiliary_model,ecm_eqnums_in_auxiliary_model,p+1) = -AR(:,:,p);
    % Write Companion matrix
    oo_.trend_component.(auxiliary_model_name).CompanionMatrix = zeros(size(B, 1)*size(B, 3));
    for i=1:p
        oo_.trend_component.(auxiliary_model_name).CompanionMatrix(1:n, (i-1)*n+(1:n)) = B(:,:,i);
        oo_.trend_component.(auxiliary_model_name).CompanionMatrix(i*n+(1:n),(i-1)*n+(1:n)) = eye(n);
    end
    oo_.trend_component.(auxiliary_model_name).CompanionMatrix(1:n, p*n+(1:n)) = B(:,:,p+1);
    M_.trend_component.(auxiliary_model_name).list_of_variables_in_companion_var = M_.endo_names(M_.trend_component.(auxiliary_model_name).lhs);
    variables_rewritten_in_levels = M_.trend_component.(auxiliary_model_name).list_of_variables_in_companion_var(ecm_eqnums_in_auxiliary_model);
    for i=1:m
        id = get_aux_variable_id(variables_rewritten_in_levels{i});
        if id
            auxinfo = M_.aux_vars(id);
            if auxinfo.type==8
                M_.trend_component.(auxiliary_model_name).list_of_variables_in_companion_var(ecm_eqnums_in_auxiliary_model(i)) = ...
                    {M_.endo_names{auxinfo.orig_index}};
            else
                error('This is a bug. Please contact the Dynare Team.')
            end
        else
            error('This is a bug. Please contact the Dynare Team.')
        end
    end
end