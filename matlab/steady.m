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

[oo_.steady_state,M_.params,info] = evaluate_steady_state(oo_.steady_state,M_,options_,oo_,~options_.steadystate.nocheck);

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


function [M,oo,info,ip,ix,ixd] = homotopy1(values, step_nbr, M, options, oo)
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
%    M              struct of model parameters
%    options        struct of options
%    oo             struct of outputs
%
% OUTPUTS
%    M              struct of model parameters
%    oo             struct of outputs
%    ip             index of parameters
%    ix             index of exogenous variables
%    ixp            index of exogenous deterministic variables

nv = size(values, 1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

if length([ip; ix; ixd]) ~= nv
    error('HOMOTOPY mode 1: incorrect variable types specified')
end

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
oldvalues = values(:,3);
ipn = find(values(:,1) == 4 & isnan(oldvalues));
oldvalues(ipn) = M.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(oldvalues));
oldvalues(ixn) = oo.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(oldvalues));
oldvalues(ixdn) = oo.exo_det_steady_state(values(ixdn, 2));

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
    old_params = M.params;
    old_exo = oo.exo_steady_state;
    old_exo_det = oo.exo_det_steady_state;
    M.params(values(ip,2)) = points(ip,i);
    oo.exo_steady_state(values(ix,2)) = points(ix,i);
    oo.exo_det_steady_state(values(ixd,2)) = points(ixd,i);

    [steady_state,M.params,info] = evaluate_steady_state(oo.steady_state,M,options,oo,~options.steadystate.nocheck);
    if info(1) == 0
        % if homotopy step is not successful, current values of steady
        % state are not modified
        oo.steady_state = steady_state;
    else
        M.params = old_params;
        oo.exo_steady_state = old_exo;
        oo.exo_det_steady_state = old_exo_det;
        break
    end
end


function homotopy2(values, step_nbr)
% Implements homotopy (mode 2) for steady-state computation.
% Only one parameter/exogenous is changed at a time.
% Computation jumps to next variable only when current variable has been
% brought to its final value.
% Variables are processed in the order in which they appear in "values".
% The problem is solved var_nbr*step_nbr times.
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

global M_ oo_ options_

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

        oo_.steady_state = evaluate_steady_state(oo_.steady_state,M_,options_,oo_,~options_.steadystate.nocheck);
    end
end


function [M,oo,info,ip,ix,ixd] = homotopy3(values, step_nbr, M, options, oo)
% Implements homotopy (mode 3) for steady-state computation.
% Tries first the most extreme values. If it fails to compute the steady
% state, the interval between initial and desired values is divided by two
% for each parameter. Every time that it is impossible to find a steady
% state, the previous interval is divided by two. When one succeed to find
% a steady state, the previous interval is multiplied by two.
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
%    step_nbr:      maximum number of steps to try before aborting
%    M              struct of model parameters
%    options        struct of options
%    oo             struct of outputs
%
% OUTPUTS
%    M              struct of model parameters
%    oo             struct of outputs
%    info           return status 0: OK, 1: failed
%    ip             index of parameters
%    ix             index of exogenous variables
%    ixp            index of exogenous deterministic variables

info = [];
tol = 1e-8;

nv = size(values,1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

if length([ip; ix; ixd]) ~= nv
    error('HOMOTOPY mode 3: incorrect variable types specified')
end

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
last_values = values(:,3);
ipn = find(values(:,1) == 4 & isnan(last_values));
last_values(ipn) = M.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(last_values));
last_values(ixn) = oo.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(last_values));
last_values(ixdn) = oo.exo_det_steady_state(values(ixdn, 2));

targetvalues = values(:,4);

%if min(abs(targetvalues - last_values)) < tol
%    error('HOMOTOPY mode 3: distance between initial and final values should be at least %e for all variables', tol)
%end
iplus = find(targetvalues > last_values);
iminus = find(targetvalues < last_values);

curvalues = last_values;
inc = (targetvalues-last_values)/2;
kplus = [];
kminus = [];
last_values = [];

disp('HOMOTOPY mode 3: launching solver at initial point...')

iter = 1;
while iter <= step_nbr

    M.params(values(ip,2)) = curvalues(ip);
    oo.exo_steady_state(values(ix,2)) = curvalues(ix);
    oo.exo_det_steady_state(values(ixd,2)) = curvalues(ixd);

    old_ss = oo.steady_state;

    [steady_state,params,info] = evaluate_steady_state(old_ss,M,options,oo,~options.steadystate.nocheck);
    if info(1) == 0
        oo.steady_state = steady_state;
        M.params = params;
        if length([kplus; kminus]) == nv
            return
        end
        if iter == 1
            disp('HOMOTOPY mode 3: successful step, now jumping to final point...')
        else
            disp('HOMOTOPY mode 3: successful step, now multiplying increment by 2...')
        end
        last_values = curvalues;
        old_params = params;
        old_exo_steady_state = oo.exo_steady_state;
        old_exo_det_steady_state = oo.exo_det_steady_state;
        inc = 2*inc;
    elseif iter == 1
        error('HOMOTOPY mode 3: can''t solve the model at 1st iteration')
    else
        disp('HOMOTOPY mode 3: failed step, now dividing increment by 2...')
        inc = inc/2;
        oo.steady_state = old_ss;
    end

    curvalues = last_values + inc;
    kplus = find(curvalues(iplus) >= targetvalues(iplus));
    curvalues(iplus(kplus)) = targetvalues(iplus(kplus));
    kminus = find(curvalues(iminus) <= targetvalues(iminus));
    curvalues(iminus(kminus)) = targetvalues(iminus(kminus));

    if max(abs(inc)) < tol
        disp('HOMOTOPY mode 3: failed, increment has become too small')
        M.params = old_params;
        oo.exo_steady_state = old_exo_steady_state;
        oo.exo_det_steady_state = old_exo_det_steady_state;
        return
    end

    iter = iter + 1;
end
disp('HOMOTOPY mode 3: failed, maximum iterations reached')
M.params = old_params;
oo.exo_steady_state = old_exo_steady_state;
oo.exo_det_steady_state = old_exo_det_steady_state;
