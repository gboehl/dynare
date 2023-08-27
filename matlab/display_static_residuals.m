function z = display_static_residuals(M_, options_, oo_,options_resid_)
% function z = display_static_residuals(M_, options_, oo_,options_resid_)
%
% Computes static residuals associated with the guess values.
%
% INPUTS
%   M:              [structure] storing the model information
%   options:        [structure] storing the options
%   oo:             [structure] storing the results
%   options_resid_:             options to resid
%
% OUTPUTS
%    z:      residuals
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2001-2023 Dynare Team
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

non_zero = (nargin > 3 && isfield(options_resid_, 'non_zero') && options_resid_.non_zero);

tags  = M_.equations_tags;
istag = 0;
if ~isempty(tags)
    istag = 1;
end

if any(imag(oo_.steady_state))
    imagrow=find(imag(oo_.steady_state));
    if ~isempty(imagrow)
        fprintf('\nresid: The initial values for the steady state of the following variables are complex:\n');
        for iter=1:length(imagrow)
            fprintf('%s\n', M_.endo_names{imagrow(iter)});
        end
    end
end

if options_.steadystate_flag
    [oo_.steady_state,M_.params] = ...
        evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,false);
end

% Compute the residuals
z = evaluate_static_model(oo_.steady_state, [oo_.exo_steady_state; ...
                                             oo_.exo_det_steady_state], ...
                          M_.params, M_, options_);

if ismember(options_.solve_algo,[10,11])
    [lb,ub,eq_index] = get_complementarity_conditions(M_,options_.ramsey_policy);
    eq_to_check=find(isfinite(lb) | isfinite(ub));
    eq_to_ignore=eq_to_check(oo_.steady_state(eq_to_check,:)<=lb(eq_to_check)+eps | oo_.steady_state(eq_to_check,:)>=ub(eq_to_check)-eps);
    z(eq_index(eq_to_ignore))=0;
    disp_string=' (accounting for MCP tags)';
else
    if istag && ~isempty(strmatch('mcp',M_.equations_tags(:,2),'exact'))
        disp_string=' (ignoring MCP tags)';
    else
        disp_string='';
    end
end

% Display the non-zero residuals if no return value
if nargout == 0
    skipline(2)
    ind = [];
    fprintf('Residuals of the static equations%s:',disp_string)
    skipline()
    any_non_zero_residual = false;
    if options_.ramsey_policy
        first_eq = M_.ramsey_orig_endo_nbr+1;
        last_eq = M_.ramsey_orig_endo_nbr+M_.ramsey_orig_eq_nbr;
    elseif options_.discretionary_policy
        first_eq = 1;
        last_eq = M_.discretionary_orig_eq_nbr;
    else
        first_eq = 1;
        last_eq = M_.orig_endo_nbr;
    end
    disp_format_tags_real=sprintf('Equation number %%%uu: %%-%us: %%14.6f\n',length(num2str(M_.eq_nbr)),size(strvcat(tags(:,3)),2)+1);
    disp_format_tags_complex=sprintf('Equation number %%%uu: %%-%us: real: %%14.6f, imaginary: %%g\n',length(num2str(M_.eq_nbr)),size(strvcat(tags(:,3)),2)+1);
    for i=first_eq:last_eq
        if abs(z(i)) < options_.solve_tolf/100
            tmp = 0;
        else
            tmp = z(i);
            any_non_zero_residual = true;
        end
        if istag
            tg = tags(cell2mat(tags(:,1)) == i,2:3); % all tags for equation i
            ind = strmatch('name', cellstr( tg(:,1) ) );
        end
        if ~(non_zero && tmp == 0)
            if ~istag || length(ind) == 0
                if ~isreal(z)
                    fprintf('Equation number %u: %g (imaginary part: %g)\n', i, real(tmp), imag(tmp));
                else
                    fprintf('Equation number %u: %g\n', i, tmp);
                end
            else
                if ~isreal(z)
                    fprintf(disp_format_tags_complex, i, tg{ind , 2}, real(tmp), imag(tmp) );
                else
                    fprintf(disp_format_tags_real, i, tg{ind , 2}, tmp);
                end
            end
        end
    end
    if non_zero && ~any_non_zero_residual
        disp('All residuals are zero')
    end
    skipline(2)
end