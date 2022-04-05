function z = resid(options_resid_)
% function z = resid(options_resid_)
%
% Computes static residuals associated with the guess values.
%
% INPUTS
%    options_resid_:   options to resid
%
% OUTPUTS
%    z:      residuals
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2022 Dynare Team
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

global M_ options_ oo_

% Properly handle the case where no input argument is given, e.g. when
% writing “z = resid;” in a .mod file (hence not using the preprocessor syntax).
non_zero = nargin > 0 && isfield(options_resid_, 'non_zero') && options_resid_.non_zero;

tags  = M_.equations_tags;
istag = 0;
if length(tags)
    istag = 1;
end

steady_state_old = oo_.steady_state;

% Keep of a copy of M_.Sigma_e
Sigma_e = M_.Sigma_e;

% Set M_.Sigma_e=0 (we evaluate the *deterministic* static model)
M_.Sigma_e = zeros(size(Sigma_e));

info = 0;
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
    [oo_.steady_state,M_.params,info] = ...
        evaluate_steady_state(oo_.steady_state,M_,options_,oo_,0);
end

% Compute the residuals
z = evaluate_static_model(oo_.steady_state, [oo_.exo_steady_state; ...
                                             oo_.exo_det_steady_state], ...
                          M_.params, M_, options_);

M_.Sigma_e = Sigma_e;


% Display the non-zero residuals if no return value
if nargout == 0
    skipline(4)
    ind = [];
    disp('Residuals of the static equations:')
    skipline()
    any_non_zero_residual = false;
    for i=1:M_.orig_eq_nbr
        if abs(z(i+M_.ramsey_eq_nbr)) < options_.solve_tolf/100
            tmp = 0;
        else
            tmp = z(i+M_.ramsey_eq_nbr);
            any_non_zero_residual = true;
        end
        if istag
            tg = tags(cell2mat(tags(:,1)) == i+M_.ramsey_eq_nbr,2:3); % all tags for equation i
            ind = strmatch('name', cellstr( tg(:,1) ) );
        end
        if ~(non_zero && tmp == 0)
            if ~istag || length(ind) == 0
                disp(['Equation number ' int2str(i) ' : ' num2str(tmp)])
            else
                t1 = tg( ind , 2 );
                s = cell2mat(t1);
                disp( ['Equation number ', int2str(i) ,' : ', num2str(tmp) ,' : ' s])
            end
        end
    end
    if non_zero && ~any_non_zero_residual
        disp('All residuals are zero')
    end
    skipline(2)
end

oo_.steady_state = steady_state_old;
