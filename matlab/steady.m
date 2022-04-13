function steady()
% function steady()
% computes and prints the steady state calculations
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2001-2019 Dynare Team
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

global M_ oo_ options_

test_for_deep_parameters_calibration(M_);

if options_.steadystate_flag && options_.homotopy_mode
    error('STEADY: Can''t use homotopy when providing a steady state external file');
end

% Keep of a copy of M_.Sigma_e
Sigma_e = M_.Sigma_e;

% Set M_.Sigma_e=0 (we compute the *deterministic* steady state)
M_.Sigma_e(:,:) = 0;

info = 0;
switch options_.homotopy_mode
  case 1
    [M_,oo_,info,ip,ix,ixd] = homotopy1(options_.homotopy_values,options_.homotopy_steps,M_,options_,oo_);
  case 2
    homotopy2(options_.homotopy_values, options_.homotopy_steps);
  case 3
    [M_,oo_,info,ip,ix,ixd] = homotopy3(options_.homotopy_values,options_.homotopy_steps,M_,options_,oo_);
end

if info(1)
    hv = options_.homotopy_values;
    skipline()
    disp('WARNING: homotopy step was not completed')
    disp('The last values for which a solution was found are:')
    for i=1:length(ip)
        fprintf('%12s %12.6f\n',char(M_.param_names(hv(ip(i),2))), ...
                M_.params(hv(ip(i),2)))
    end
    for i=1:length(ix)
        fprintf('%12s %12.6f\n',char(M_.exo_names(hv(ix(i),2))), ...
                oo_.exo_steady_state(hv(ix(i),2)))
    end
    for i=1:length(ixd)
        fprintf('%12s %12.6f\n',char(M_.exo_det_names(hv(ixd(i),2))), ...
                oo_.exo_det_steady_state(hv(ixd(i),2)))
    end

    if options_.homotopy_force_continue
        disp('Option homotopy_continue is set, so I continue ...')
    else
        error('Homotopy step failed')
    end
end

[oo_.steady_state,M_.params,info] = steady_(M_,options_,oo_);

if info(1) == 0
    if ~options_.noprint
        disp_steady_state(M_,oo_,options_);
    end
else
    if ~options_.noprint
        if ~isempty(oo_.steady_state)
            resid;
        else
            skipline()
            disp('Residuals of the static equations cannot be computed because the steady state routine returned an empty vector.')
            skipline()
        end
    end
    if options_.debug
        fprintf('\nThe steady state computation failed. It terminated with the following values:\n')
        for i=1:M_.orig_endo_nbr
            fprintf('%s \t\t %g\n', M_.endo_names{i}, oo_.steady_state(i));
        end
    end
    print_info(info,options_.noprint, options_);
end

M_.Sigma_e = Sigma_e;
