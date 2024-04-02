function [expression, growthneutralitycorrection] = write_pac_mce_expectations(eqname, expectationmodelname, auxname)

% Prints the expansion of the PAC_EXPECTATION term in files.
%
% INPUTS
% - eqname                         [char]       Name of the equation.
% - epxpectationmodelname          [char]       Name of the expectation model.
% - iscrlf                         [logical]    Adds carriage return after each additive term if true (default is false).
%
% OUTPUTS
% - expression                     [char]       Unrolled expectation expression.
% - growthneutralitycorrection     [char]

% Copyright © 2019-2024 Dynare Team
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

expectationmodel = M_.pac.(expectationmodelname);

targetid = expectationmodel.ec.vars((expectationmodel.ec.istarget==true));
alphaid = expectationmodel.mce.alpha;
betaid = expectationmodel.discount_index;

target = M_.endo_names{targetid};
transformations = {};

if isauxiliary(target)
    ida = get_aux_variable_id(target);
    op = 0;
    while ida
        op = op+1;
        if isequal(M_.aux_vars(ida).type, 8)
            transformations(op) = {'diff'};
            target = M_.endo_names{M_.aux_vars(ida).orig_index};
            ida = get_aux_variable_id(target);
        elseif isequal(M_.aux_vars(ida).type, 10)
            transformations(op) = {M_.aux_vars(ida).unary_op};
            target = M_.endo_names{M_.aux_vars(ida).orig_index};
            ida = get_aux_variable_id(target);
        else
            error('This case is not implemented.')
        end
    end
end


% In Brayton, Davis and Tulip (2000) the formula for the PAC expectation under Model Consistent Expectations is given as:
%
%       ᵐ                              ᵐ⁻¹ ᵐ⁻¹
%  Zₜ =  ∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ + A(1) [ 𝛥 yₜ − ∑   ∑   𝛼ⱼ₊₁𝛽ʲ⁺¹𝛥 yₜ₊ₖ ]
%       ᵢ₌₁                            ₖ₌₁ ⱼ₌ₖ
%
% where yₜ is the target in the PAC equation, see equation 10 or A.92 in the appendix. The sign before the sum of expected
% Z is incorrect, in the follwing code we define Z as:
%
%        ᵐ                              ᵐ⁻¹ ᵐ⁻¹
%  Zₜ =  -∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ + A(1) [ 𝛥 yₜ − ∑   ∑   𝛼ⱼ₊₁𝛽ʲ⁺¹𝛥 yₜ₊ₖ ]
%        ᵢ₌₁                            ₖ₌₁ ⱼ₌ₖ

expression = '';
A1 = '1';

% First loop to build
%
%     ᵐ
%    -∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ
%    ᵢ₌₁
%
%              ᵐ
% and A(1)=1 + ∑ αᵢ
%             ᵢ₌₁

for i=1:length(alphaid)
    expression = sprintf('%s-%s*(%s^%i)*%s(%i)', expression, ...
                         M_.param_names{alphaid(i)}, ...
                         M_.param_names{betaid}, ...
                         i, auxname, i);
    A1 = sprintf('%s+%s', A1, M_.param_names{alphaid(i)});
end

% Write
%
%     ᵐ
%    -∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ + A(1)
%    ᵢ₌₁


expression = sprintf('%s+(%s)', expression, A1);

% Write
%
%     ᵐ
%    -∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ + A(1) [Δyₜ
%    ᵢ₌₁

if isempty(transformations)
    expression = sprintf('%s*(diff(%s)', expression, target);
else
    % Typically the target may be the log of a variable or the first difference of the log of a variable.
    variable = target;
    for k=length(transformations):-1:1
        variable = sprintf('%s(%s)', transformations{k}, variable);
    end
    expression = sprintf('%s*(diff(%s)', expression, variable);
end

% Second loop (for the double sum) to write
%
%        ᵐ                              ᵐ⁻¹ ᵐ⁻¹
%  Zₜ =  -∑  𝛼ᵢ 𝛽ⁱ⁺¹ Zₜ₊ᵢ + A(1) [ 𝛥 yₜ − ∑   ∑   𝛼ⱼ₊₁𝛽ʲ⁺¹𝛥 yₜ₊ₖ
%        ᵢ₌₁                            ₖ₌₁ ⱼ₌ₖ

for i=1:length(alphaid)-1
    tmp = sprintf('%s*%s^%i', M_.param_names{alphaid(i+1)}, M_.param_names{betaid}, i+1);
    for j=i+1:length(alphaid)-1
        tmp = sprintf('%s+%s*%s^%i', tmp, M_.param_names{alphaid(j+1)}, M_.param_names{betaid}, j+1);
    end
    if isempty(transformations)
        expression = sprintf('%s-(%s)*diff(%s(%i))', expression, tmp, target, i);
    else
        variable = sprintf('%s(%i)', target, i);
        for k=length(transformations):-1:1
            variable = sprintf('%s(%s)', transformations{k}, variable);
        end
        expression = sprintf('%s-(%s)*diff(%s)', expression, tmp, variable);
    end
end

% Close brackets.
expression = sprintf('%s)', expression);

% Add growth neutrality correction if required.
if isfield(expectationmodel, 'growth_neutrality_param_index')
    if numel(expectationmodel.growth_linear_comb) == 1
        growthneutralitycorrection = sprintf('%s*%s', M_.param_names{expectationmodel.growth_neutrality_param_index}, expectationmodel.growth_str);
    else
        growthneutralitycorrection = sprintf('%s*(%s)', M_.param_names{expectationmodel.growth_neutrality_param_index}, expectationmodel.growth_str);
    end
else
    growthneutralitycorrection = '';
end
