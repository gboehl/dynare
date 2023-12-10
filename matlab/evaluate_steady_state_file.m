function [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M_,options_,steady_state_checkflag)
% function [ys,params1,info] = evaluate_steady_state_file(ys_init,exo_ss,M_,options_,steady_state_checkflag)
% Evaluates steady state files
%
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state (incl. deterministic exogenous)
%   M_                        struct           model parameters
%   options_                  struct           options
%   steady_state_checkflag    boolean          indicator whether to check steady state returned
% OUTPUTS
%   ys                        vector           steady state
%   params1                   vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

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

params = M_.params;
info = 0;

fname = M_.fname;

if options_.steadystate_flag == 1
    % old format
    h_steadystate = str2func([fname '_steadystate']);
    [ys,params1,check] = h_steadystate(ys_init, exo_ss,M_,options_);
else % steadystate_flag == 2
     % new format
    h_steadystate = str2func([fname '.steadystate']);
    [ys,params1,check] = h_steadystate(ys_init, exo_ss, params);
end

if check
    info(1) = 19;
    info(2) = NaN;
    return
end

if M_.param_nbr > 0
    updated_params_flag = max(abs(params1-params)) > 1e-12 ...
        || ~isequal(isnan(params1),isnan(params)); %checks whether numbers or NaN changed
else
    updated_params_flag = 0;
end

h_set_auxiliary_variables = str2func([M_.fname '.set_auxiliary_variables']);

if  isnan(updated_params_flag) || (updated_params_flag  && any(isnan(params(~isnan(params))-params1(~isnan(params))))) %checks if new NaNs were added
    info(1) = 24;
    info(2) = NaN;
    if M_.set_auxiliary_variables
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
    end
    return
end

if updated_params_flag && ~isreal(params1)
    info(1) = 23;
    if ~isoctave
        info(2) = sum(imag(params).^2,'omitnan');
    else
        info(2) = nansum(imag(params).^2);
    end
    if M_.set_auxiliary_variables
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
    end
    return
end

if updated_params_flag
    params = params1;
end

% adding values for auxiliary variables
if ~isempty(M_.aux_vars) && ~options_.ramsey_policy
    if M_.set_auxiliary_variables
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
    end
end

if steady_state_checkflag
    % Check whether the steady state obtained from the _steadystate file is a steady state.
    [residuals, check] = evaluate_static_model(ys, exo_ss, params, M_, options_);
    if check
        info(1) = 19;
        info(2) = check; % to be improved
        return
    end
    if max(abs(residuals)) > options_.solve_tolf
        info(1) = 19;
        info(2) = residuals'*residuals;
        return
    end
    if any(isnan(residuals))
        info(1) = 22;
        return
    end
end
