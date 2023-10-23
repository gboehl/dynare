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

global M_ oo_ options_

test_for_deep_parameters_calibration(M_);

% Keep of a copy of M_.Sigma_e
Sigma_e = M_.Sigma_e;

% Set M_.Sigma_e=0 (we compute the *deterministic* steady state)
M_.Sigma_e(:,:) = 0;

if ~ismember(options_.homotopy_mode, [0 1 2 3])
    error('STEADY: invalid value for homotopy_mode option')
end

if isfield(options_, 'homotopy_values') && options_.homotopy_mode == 0
    warning('STEADY: a homotopy_setup block is present but homotopy will not be performed because homotopy_mode option is equal to 0')
end

if options_.homotopy_mode ~= 0
    if ~isfield(options_, 'homotopy_values')
        error('STEADY: a homotopy_setup block must be present when the homotopy_mode option is specified')
    end

    if options_.steadystate_flag
        error('STEADY: Can''t use homotopy when providing a steady state external file');
    end

    hv = options_.homotopy_values;

    if any(hv(:,1)~=1 & hv(:,1)~=2 & hv(:,1)~=4)
        % Already checked by the preprocessor, but let’s stay on the safe side
        error('HOMOTOPY_SETUP: incorrect variable types specified')
    end

    % If the “from_initval_to_endval” option was passed to the “homotopy_setup” block, add the relevant homotopy information
    if options_.homotopy_from_initval_to_endval
        if isempty(oo_.initial_exo_steady_state)
            error('HOMOTOPY_SETUP: the from_initval_to_endval option cannot be used without an endval block')
        end
        for i = 1:M_.exo_nbr
            if ~any(hv(:,1)==1 & hv(:,2)==i) % Do not overwrite information manually specified by the user
                hv = vertcat(hv, [ 1 i oo_.initial_exo_steady_state(i) oo_.exo_steady_state(i)]);
            end
        end
    end

    homotopy_func = str2func(['homotopy' num2str(options_.homotopy_mode)]);
    [M_,oo_,errorcode] = homotopy_func(hv, options_.homotopy_steps, M_, options_, oo_);

    if errorcode
        if errorcode == 2
            disp('WARNING: homotopy failed at the first iteration (for the starting values)')
        else % errorcode == 1: print last successful point
            ip = find(hv(:,1) == 4); % Parameters
            ix = find(hv(:,1) == 1); % Exogenous
            ixd = find(hv(:,1) == 2); % Exogenous deterministic
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
        end

        if options_.homotopy_force_continue
            disp('Option homotopy_continue is set, so I continue ...')
        else
            error('Homotopy step failed')
        end
    end
end

[oo_.steady_state,M_.params,info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);

if info(1) == 0
    if ~options_.noprint
        disp_steady_state(M_,oo_,options_);
    end
else
    if ~options_.noprint
        if ~isempty(oo_.steady_state)
            display_static_residuals(M_, options_, oo_);
        else
            skipline()
            disp('Residuals of the static equations cannot be computed because the steady state routine returned an empty vector.')
            skipline()
        end
    end
    if options_.debug
        fprintf('\nsteady: The steady state computation failed. It terminated with the following values:\n')
        if ~isreal(oo_.steady_state)
            format_string=sprintf('%%-%us= %%g%%+gi\n',size(strvcat(M_.endo_names),2)+1);
        else
            format_string=sprintf('%%-%us= %%14.6f\n',size(strvcat(M_.endo_names),2)+1);
        end
        for i=1:M_.orig_endo_nbr
            if ~isreal(oo_.steady_state)
                fprintf(format_string, M_.endo_names{i}, real(oo_.steady_state(i)),imag(oo_.steady_state(i)));
            else
                fprintf(format_string, M_.endo_names{i}, oo_.steady_state(i));
            end
        end
    end
    print_info(info,options_.noprint, options_);
end

M_.Sigma_e = Sigma_e;


function [M_,oo_,errorcode] = homotopy1(values, step_nbr, M_, options_, oo_)
% Implements homotopy (mode 1) for steady-state computation.
% The multi-dimensional vector going from the set of initial values
% to the set of final values is divided in as many sub-vectors as
% there are steps, and the problem is solved as many times.
%
% INPUTS
%    values:        a matrix with 4 columns, representing the content of
%                   homotopy_setup block, with one variable per line.
%                   Column 1 is variable type (1 for exogenous, 2 for
%                   exogenous deterministic, 4 for parameters)
%                   Column 2 is symbol integer identifier.
%                   Column 3 is initial value, and column 4 is final value.
%                   Column 3 can contain NaNs, in which case previous
%                   initialization of variable will be used as initial value.
%    step_nbr:      number of steps for homotopy
%    M_             struct of model parameters
%    options_       struct of options
%    oo_            struct of outputs
%
% OUTPUTS
%    M_             struct of model parameters
%    oo_            struct of outputs
%    errorcode      0 in case of success
%                   1 if some homotopy steps were successful but it was not
%                     possible to go up to 100%; in that case, parameters in
%                     M_.params and exogenous in oo_ are left to the last
%                     successful point
%                   2 if it wasn’t possible to compute a solution for the
%                     starting values

nv = size(values, 1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
oldvalues = values(:,3);
ipn = find(values(:,1) == 4 & isnan(oldvalues));
oldvalues(ipn) = M_.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(oldvalues));
oldvalues(ixn) = oo_.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(oldvalues));
oldvalues(ixdn) = oo_.exo_det_steady_state(values(ixdn, 2));

points = zeros(nv, step_nbr+1);
for i = 1:nv
    if (oldvalues(i) ~= values(i, 4))
        points(i,:) = oldvalues(i):(values(i,4)-oldvalues(i))/step_nbr:values(i,4);
    else
        points(i,:) = values(i,4);
    end
end

for i=1:step_nbr+1
    disp([ 'HOMOTOPY mode 1: computing step ' int2str(i-1) '/' int2str(step_nbr) '...' ])
    M_.params(values(ip,2)) = points(ip,i);
    oo_.exo_steady_state(values(ix,2)) = points(ix,i);
    oo_.exo_det_steady_state(values(ixd,2)) = points(ixd,i);
    [oo_.steady_state,M_.params,info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);
    if info(1)
        if i == 1
            errorcode = 2;
        else
            M_.params = last_successful_params;
            oo_.exo_steady_state = last_successful_exo;
            oo_.exo_det_steady_state = last_successful_exo_det;
            errorcode = 1;
        end
        return
    end
    last_successful_params = M_.params;
    last_successful_exo = oo_.exo_steady_state;
    last_successful_exo_det = oo_.exo_det_steady_state;
end

errorcode = 0;


function [M_, oo_, errorcode] = homotopy2(values, step_nbr, M_, options_, oo_)
% Implements homotopy (mode 2) for steady-state computation.
% Only one parameter/exogenous is changed at a time.
% Computation jumps to next variable only when current variable has been
% brought to its final value.
% Variables are processed in the order in which they appear in "values".
% The problem is solved var_nbr*step_nbr times.
%
% See homotopy1 for the description of inputs and outputs.

nv = size(values, 1);

oldvalues = values(:,3);

% Initialize all variables with initial value, or the reverse...
for i = 1:nv
    switch values(i,1)
      case 1
        if isnan(oldvalues(i))
            oldvalues(i) = oo_.exo_steady_state(values(i,2));
        else
            oo_.exo_steady_state(values(i,2)) = oldvalues(i);
        end
      case 2
        if isnan(oldvalues(i))
            oldvalues(i) = oo_.exo_det_steady_state(values(i,2));
        else
            oo_.exo_det_steady_state(values(i,2)) = oldvalues(i);
        end
      case 4
        if isnan(oldvalues(i))
            oldvalues(i) = M_.params(values(i,2));
        else
            M_.params(values(i,2)) = oldvalues(i);
        end
      otherwise
        error('HOMOTOPY mode 2: incorrect variable types specified')
    end
end

if any(oldvalues == values(:,4))
    error('HOMOTOPY mode 2: initial and final values should be different')
end

% Actually do the homotopy
for i = 1:nv
    switch values(i,1)
      case 1
        varname = M_.exo_names{values(i,2)};
      case 2
        varname = M_.exo_det_names{values(i,2)};
      case 4
        varname = M_.param_names{values(i,2)};
    end
    for v = oldvalues(i):(values(i,4)-oldvalues(i))/step_nbr:values(i,4)
        switch values(i,1)
          case 1
            oo_.exo_steady_state(values(i,2)) = v;
          case 2
            oo_.exo_det_steady_state(values(i,2)) = v;
          case 4
            M_.params(values(i,2)) = v;
        end

        disp([ 'HOMOTOPY mode 2: lauching solver with ' varname ' = ' num2str(v) ' ...'])

        [oo_.steady_state, M_.params, info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);
        if info(1)
            if i == 1 && v == oldvalues(1)
                errorcode = 2;
            else
                M_.params = last_successful_params;
                oo_.exo_steady_state = last_successful_exo;
                oo_.exo_det_steady_state = last_successful_exo_det;
                errorcode = 1;
            end
            return
        end
        last_successful_params = M_.params;
        last_successful_exo = oo_.exo_steady_state;
        last_successful_exo_det = oo_.exo_det_steady_state;
    end
end

errorcode = 0;


function [M_,oo_,errorcode] = homotopy3(values, step_nbr, M_, options_, oo_)
% Implements homotopy (mode 3) for steady-state computation.
% Tries first the most extreme values. If it fails to compute the steady
% state, the interval between initial and desired values is divided by two
% for each parameter. Every time that it is impossible to find a steady
% state, the previous interval is divided by two. When one succeed to find
% a steady state, the previous interval is multiplied by two.
%
% See homotopy1 for the description of inputs and outputs.

info = [];
tol = 1e-8;

nv = size(values,1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
oldvalues = values(:,3);
ipn = find(values(:,1) == 4 & isnan(oldvalues));
oldvalues(ipn) = M_.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(oldvalues));
oldvalues(ixn) = oo_.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(oldvalues));
oldvalues(ixdn) = oo_.exo_det_steady_state(values(ixdn, 2));

targetvalues = values(:,4);

iplus = find(targetvalues > oldvalues);
iminus = find(targetvalues < oldvalues);

curvalues = oldvalues;
inc = (targetvalues-oldvalues)/2;
kplus = [];
kminus = [];
last_successful_values = [];

disp('HOMOTOPY mode 3: launching solver at initial point...')

iter = 1;
while iter <= step_nbr

    M_.params(values(ip,2)) = curvalues(ip);
    oo_.exo_steady_state(values(ix,2)) = curvalues(ix);
    oo_.exo_det_steady_state(values(ixd,2)) = curvalues(ixd);

    old_ss = oo_.steady_state;

    [steady_state,params,info] = evaluate_steady_state(old_ss,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);
    if info(1) == 0
        oo_.steady_state = steady_state;
        M_.params = params;
        if length([kplus; kminus]) == nv
            errorcode = 0;
            return
        end
        if iter == 1
            disp('HOMOTOPY mode 3: successful step, now jumping to final point...')
        else
            disp('HOMOTOPY mode 3: successful step, now multiplying increment by 2...')
        end
        last_successful_values = curvalues;
        last_successful_params = params;
        last_successful_exo_steady_state = oo_.exo_steady_state;
        last_successful_exo_det_steady_state = oo_.exo_det_steady_state;
        inc = 2*inc;
    elseif iter == 1
        errorcode = 2;
        return
    else
        disp('HOMOTOPY mode 3: failed step, now dividing increment by 2...')
        inc = inc/2;
        oo_.steady_state = old_ss;
    end

    curvalues = last_successful_values + inc;
    kplus = find(curvalues(iplus) >= targetvalues(iplus));
    curvalues(iplus(kplus)) = targetvalues(iplus(kplus));
    kminus = find(curvalues(iminus) <= targetvalues(iminus));
    curvalues(iminus(kminus)) = targetvalues(iminus(kminus));

    if max(abs(inc)) < tol
        disp('HOMOTOPY mode 3: failed, increment has become too small')
        M_.params = last_successful_params;
        oo_.exo_steady_state = last_successful_exo_steady_state;
        oo_.exo_det_steady_state = last_successful_exo_det_steady_state;
        errorcode = 1;
        return
    end

    iter = iter + 1;
end
disp('HOMOTOPY mode 3: failed, maximum iterations reached; you may want to increase the homotopy_steps option')
M_.params = last_successful_params;
oo_.exo_steady_state = last_successful_exo_steady_state;
oo_.exo_det_steady_state = last_successful_exo_det_steady_state;
errorcode = 1;
