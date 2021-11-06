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
number_of_pac_eq = size(pacmodel.tag_map, 1);
if number_of_pac_eq==1
    fprintf('PAC model %s is used in equation %s.\n', pacname, pacmodel.tag_map{1,1});
else
    fprintf('PAC model %s is used in %u equation(s):\n', pacname, number_of_pac_eq);
    skipline()
    for i=1:number_of_pac_eq
        fprintf('    -  %s\n', pacmodel.tag_map{i,1});
    end
end
skipline()

equations = pacmodel.equations;

for e=1:number_of_pac_eq
    % Get PAC equation tag
    eqtag = pacmodel.tag_map{e,2};
    % Get PAC equation
    pac_equation = equations.(eqtag);
    % Get Error correction and autoregressive parameters in PAC equation
    params = NaN(2+pac_equation.max_lag, 1);
    params(1) = M_.params(pac_equation.ec.params);
    params(1+(1:pac_equation.max_lag)) = M_.params(pac_equation.ar.params);
    params(end) = M_.params(pacmodel.discount_index);
    [G, alpha, beta] = buildGmatrixWithAlphaAndBeta(params);
    M_.params(pac_equation.mce.alpha) = flip(alpha);
    if isfield(pacmodel, 'growth_neutrality_param_index')
        if isfield(equations.(eqtag), 'non_optimizing_behaviour')
            gamma = M_.params(equations.(eqtag).share_of_optimizing_agents_index);
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
        if isfield(equations.(eqtag), 'optim_additive')
            % Exogenous variables are present in the λ part (optimizing agents).
            tmp0 = 0;
            for i=1:length(equations.(eqtag).optim_additive.params)
                if isnan(equations.(eqtag).optim_additive.params(i)) && islogical(equations.(eqtag).optim_additive.bgp{i}) && equations.(eqtag).optim_additive.bgp{i}
                    tmp0 = tmp0 + equations.(eqtag).optim_additive.scaling_factor(i);
                elseif ~isnan(equations.(eqtag).optim_additive.params(i)) && islogical(equations.(eqtag).optim_additive.bgp{i}) && equations.(eqtag).optim_additive.bgp{i}
                    tmp0 = tmp0 + M_.params(equations.(eqtag).optim_additive.params(i))*equations.(eqtag).optim_additive.scaling_factor(i);
                elseif ~islogical(equations.(eqtag).optim_additive.bgp{i})
                    error('It is not possible to provide a value for the mean of an exogenous variable appearing in the optimal part of the PAC equation.')
                end
            end
            cc = cc - tmp0;
        end
        if gamma<1
            if isfield(equations.(eqtag), 'non_optimizing_behaviour') && isfield(equations.(eqtag).non_optimizing_behaviour, 'params')
                % Exogenous variables are present in the 1-λ part (rule of thumb agents).
                tmp0 = 0;
                tmp1 = 0;
                for i=1:length(equations.(eqtag).non_optimizing_behaviour.params)
                    if isnan(equations.(eqtag).non_optimizing_behaviour.params(i)) && islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && equations.(eqtag).non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i);
                    elseif ~isnan(equations.(eqtag).non_optimizing_behaviour.params(i)) && islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && equations.(eqtag).non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + M_.params(equations.(eqtag).non_optimizing_behaviour.params(i))*equations.(eqtag).non_optimizing_behaviour.scaling_factor(i);
                    elseif ~islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnumeric(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnan(equations.(eqtag).non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i)*equations.(eqtag).non_optimizing_behaviour.bgp{i};
                    elseif ~islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnumeric(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && ~isnan(equations.(eqtag).non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i)*equations.(eqtag).non_optimizing_behaviour.params(i)*equations.(eqtag).non_optimizing_behaviour.bgp{i};
                    end
                end
                cc = cc - (1-gamma)*tmp0/gamma;
                ll = -(1.0-gamma)*tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
            end
        end
        if isfield(equations.(eqtag), 'additive')
            % Exogenous variables are present outside of the λ and (1-λ) parts (or we have exogenous variables in a "pure" PAC equation.
            tmp0 = 0;
            tmp1 = 0;
            for i=1:length(equations.(eqtag).additive.params)
                if isnan(equations.(eqtag).additive.params(i)) && islogical(equations.(eqtag).additive.bgp{i}) && equations.(eqtag).additive.bgp{i}
                    tmp0 = tmp0 + equations.(eqtag).additive.scaling_factor(i);
                elseif ~isnan(equations.(eqtag).additive.params(i)) && islogical(equations.(eqtag).additive.bgp{i}) && equations.(eqtag).additive.bgp{i}
                    tmp0 = tmp0 + M_.params(equations.(eqtag).additive.params(i))*equations.(eqtag).additive.scaling_factor(i);
                elseif ~islogical(equations.(eqtag).additive.bgp{i}) && isnumeric(equations.(eqtag).additive.bgp{i}) && isnan(equations.(eqtag).additive.params(i))
                    tmp1 = tmp1 + equations.(eqtag).additive.scaling_factor(i)*equations.(eqtag).additive.bgp{i};
                elseif ~islogical(equations.(eqtag).additive.bgp{i}) && isnumeric(equations.(eqtag).additive.bgp{i}) && ~isnan(equations.(eqtag).additive.params(i))
                    tmp1 = tmp1 + equations.(eqtag).additive.scaling_factor(i)*equations.(eqtag).additive.params(i)*equations.(eqtag).additive.bgp{i};
                end
            end
            cc = cc - tmp0/gamma;
            ll = ll - tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
        end
        M_.params(pacmodel.growth_neutrality_param_index) = cc; % Multiplies the variable or expression provided though the growth option in command pac_model.
    end
end