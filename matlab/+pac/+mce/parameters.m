function parameters(pacname)

% Updates the parameters of a PAC Model Consistent Expectations.
%
% INPUTS
% - pacname       [string]    Name of the pac equation.
%
% OUTPUTS
% - none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2019-2021 Dynare Team
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

global M_

% Check that the first input is a row character array.
if ~isrow(pacname)==1 || ~ischar(pacname)
    error('Input argument must be a row character array!')
end

% Check the name of the PAC model.
if ~isfield(M_.pac, pacname)
    error('PAC model %s is not defined in the model block!', pacname)
end

% Get PAC model description
pacmodel = M_.pac.(pacname);

% Check that we are dealing with PAC/MCE
if ~pacmodel.model_consistent_expectations
    error('This function is to be used only for PAC Model Consistent Expectations.')
end

% Show the equations where this PAC model is used.
fprintf('PAC model %s is used in equation %s.\n', pacname, pacmodel.eq_name);
skipline()

% Get Error correction and autoregressive parameters in PAC equation
params = NaN(2+pacmodel.max_lag, 1);
params(1) = M_.params(pacmodel.ec.params);
params(1+(1:pacmodel.max_lag)) = M_.params(pacmodel.ar.params);
params(end) = M_.params(pacmodel.discount_index);
[G, alpha, beta] = buildGmatrixWithAlphaAndBeta(params);
M_.params(pacmodel.mce.alpha) = flip(alpha);
if isfield(pacmodel, 'growth_neutrality_param_index')
    if isfield(pacmodel, 'non_optimizing_behaviour')
        gamma = M_.params(pacmodel.share_of_optimizing_agents_index);
    else
        gamma = 1.0;
    end
    A = [alpha; 1];
    A_1 = polyval(A, 1.0);
    A_b = polyval(A, beta);
    m = length(alpha);
    d = A_1*A_b*(iota(m, m)'*inv((eye(m)-G)*(eye(m)-G))*iota(m, m));
    cc = 1/gamma-(sum(params(2:end-1))+d);
    ll = 0;
    if isfield(pacmodel, 'optim_additive')
        % Exogenous variables are present in the λ part (optimizing agents).
        tmp0 = 0;
        for i=1:length(pacmodel.optim_additive.params)
            if isnan(pacmodel.optim_additive.params(i)) && islogical(pacmodel.optim_additive.bgp{i}) && pacmodel.optim_additive.bgp{i}
                tmp0 = tmp0 + pacmodel.optim_additive.scaling_factor(i);
            elseif ~isnan(pacmodel.optim_additive.params(i)) && islogical(pacmodel.optim_additive.bgp{i}) && pacmodel.optim_additive.bgp{i}
                tmp0 = tmp0 + M_.params(pacmodel.optim_additive.params(i))*equations.(eqtag).optim_additive.scaling_factor(i);
            elseif ~islogical(pacmodel.optim_additive.bgp{i})
                error('It is not possible to provide a value for the mean of an exogenous variable appearing in the optimal part of the PAC equation.')
            end
        end
        cc = cc - tmp0;
    end
    if gamma<1
        if isfield(pacmodel, 'non_optimizing_behaviour') && isfield(pacmodel.non_optimizing_behaviour, 'params')
            % Exogenous variables are present in the 1-λ part (rule of thumb agents).
            tmp0 = 0;
            tmp1 = 0;
            for i=1:length(pacmodel.non_optimizing_behaviour.params)
                if isnan(pacmodel.non_optimizing_behaviour.params(i)) && islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && pacmodel.non_optimizing_behaviour.bgp{i}
                    tmp0 = tmp0 + pacmodel.non_optimizing_behaviour.scaling_factor(i);
                elseif ~isnan(pacmodel.non_optimizing_behaviour.params(i)) && islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && pacmodel.non_optimizing_behaviour.bgp{i}
                    tmp0 = tmp0 + M_.params(pacmodel.non_optimizing_behaviour.params(i))*pacmodel.non_optimizing_behaviour.scaling_factor(i);
                elseif ~islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && isnumeric(pacmodel.non_optimizing_behaviour.bgp{i}) && isnan(pacmodel.non_optimizing_behaviour.params(i))
                    tmp1 = tmp1 + pacmodel.non_optimizing_behaviour.scaling_factor(i)*pacmodel.non_optimizing_behaviour.bgp{i};
                elseif ~islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && isnumeric(pacmodel.non_optimizing_behaviour.bgp{i}) && ~isnan(pacmodel.non_optimizing_behaviour.params(i))
                    tmp1 = tmp1 + pacmodel.non_optimizing_behaviour.scaling_factor(i)*pacmodel.non_optimizing_behaviour.params(i)*pacmodel.non_optimizing_behaviour.bgp{i};
                end
            end
            cc = cc - (1-gamma)*tmp0/gamma;
            ll = -(1.0-gamma)*tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
        end
    end
    if isfield(pacmodel, 'additive')
        % Exogenous variables are present outside of the λ and (1-λ) parts (or we have exogenous variables in a "pure" PAC equation.
        tmp0 = 0;
        tmp1 = 0;
        for i=1:length(pacmodel.additive.params)
            if isnan(pacmodel.additive.params(i)) && islogical(pacmodel.additive.bgp{i}) && pacmodel.additive.bgp{i}
                tmp0 = tmp0 + pacmodel.additive.scaling_factor(i);
            elseif ~isnan(pacmodel.additive.params(i)) && islogical(pacmodel.additive.bgp{i}) && pacmodel.additive.bgp{i}
                tmp0 = tmp0 + M_.params(pacmodel.additive.params(i))*pacmodel.additive.scaling_factor(i);
            elseif ~islogical(pacmodel.additive.bgp{i}) && isnumeric(pacmodel.additive.bgp{i}) && isnan(pacmodel.additive.params(i))
                tmp1 = tmp1 + pacmodel.additive.scaling_factor(i)*pacmodel.additive.bgp{i};
            elseif ~islogical(pacmodel.additive.bgp{i}) && isnumeric(pacmodel.additive.bgp{i}) && ~isnan(pacmodel.additive.params(i))
                tmp1 = tmp1 + pacmodel.additive.scaling_factor(i)*pacmodel.additive.params(i)*pacmodel.additive.bgp{i};
            end
        end
        cc = cc - tmp0/gamma;
        ll = ll - tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
    end
    M_.params(pacmodel.growth_neutrality_param_index) = cc; % Multiplies the variable or expression provided though the growth option in command pac_model.
end