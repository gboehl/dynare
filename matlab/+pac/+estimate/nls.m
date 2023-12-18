function nls(eqname, params, data, range, optimizer, varargin)

% Estimates the parameters of a PAC equation by Nonlinear Least Squares.
%
% INPUTS
% - eqname       [string]    Name of the pac equation.
% - params       [struct]    Describes the parameters to be estimated.
% - data         [dseries]   Database for the estimation
% - range        [dates]     Range of dates for the estimation.
% - optimizer    [string]    Set the optimization algorithm. Possible values
%                            are 'fmincon', 'fminunc', 'csminwel',
%                            'fminsearch', 'simplex' and
%                            'annealing'. Default is 'csminwel'.
% - varargin                 List of key/value pairs, each key must be a
%                            string and values can be strings or real numbers.
%
% OUTPUTS
% - none
%
% REMARKS
% [1] The estimation results are printed in the command line window, and the
%     parameters are updated accordingly in M_.params.
% [2] The second input is a structure. Each fieldname corresponds to the
%     name of an estimated parameter, the value of the field is the initial
%     guess used for the estimation (by NLS).
% [3] The third input is a dseries object which must at least contain all
%     the variables appearing in the estimated equation. The residual of the
%     equation must have NaN values in the object.
% [4] It is assumed that the residual is additive.
% [5] If the used optimization routines handles inequalities over the
%     estimated parameters, we impose the positivity of the error correction parameter.
%
% EXAMPLE
%
% pac.estimate.nls('zpac', eparams, dbase, dates('2003Q1'):dates('2016Q3'),'fminunc', 'MaxIter', 50);
%
% where zpac is the name of the PAC equation, eparams is a structure
% containing the guess values for the estimated parameters (in each field),
% dbase is a dseries object containing the data,
% dates('2003Q1'):dates('2016Q3') is the range of the sample used for
% estimation, 'fminunc' is the name of the optimization algorithm (this one
% is available only if the matylab optimization toolbox is installed), the
% remaining inputs are the options (key/value) passed to the optimizers.

% Copyright © 2018-2023 Dynare Team
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

is_gauss_newton = false;
is_lsqnonlin = false;
if nargin>4 && (isequal(optimizer, 'GaussNewton') || isequal(optimizer, 'lsqnonlin'))
    switch optimizer
      case 'GaussNewton'
        is_gauss_newton = true;
      case 'lsqnonlin'
        is_lsqnonlin = true;
    end
end

[pacmodl, lhs, rhs, pnames, enames, xnames, ~, pid, eid, xid, ~, ipnames_, params, data, islaggedvariables] = ...
    pac.estimate.init(M_, oo_, eqname, params, data, range);

dLHS = M_.aux_vars(strmatch(lhs,M_.endo_names, 'exact')==[M_.aux_vars(:).endo_index]).orig_expr;

% Check that the error correction term is correct.
if M_.pac.(pacmodl).ec.istarget(2)
    error(['\nThe error correction term in PAC equation (%s) is not correct.\nThe ' ...
           'error correction term should be the difference between a trend\n' ...
           'and the level of the endogenous variable.'], pacmodl);
end

%
% Rewrite and print the equation.
%

[rhs, lhs] = rewrite_equation_with_tables(rhs, lhs, islaggedvariables, pnames, enames, xnames, pid, eid, xid, data);

% Create a routine for evaluating the residuals of the nonlinear model
write_residuals_routine(lhs, rhs, eqname, ipnames_, M_, pacmodl);

% Create a routine for evaluating the sum of squared residuals of the nonlinear model
write_ssr_routine(lhs, rhs, eqname, ipnames_, M_, pacmodl);

% Copy (sub)sample data in a matrix.
DATA = data([range(1)-1, range]).data;

% Create a function handle returning the sum of square residuals for a given vector of parameters.
ssrfun = @(p) feval([M_.fname '.ssr_' eqname], p, DATA, M_, oo_);

% Create a function handle returning the sum of square residuals for a given vector of parameters.
resfun = @(p) feval([M_.fname '.r_' eqname], p, DATA, M_, oo_);

% Set initial condition.
params0 = cell2mat(struct2cell(params));

% Set defaults for optimizers
bounds = [];
parameter_names = [];

% Set optimizer routine.
if nargin<5 || isempty(optimizer)
    % Use csminwel by default.
    minalgo = 4;
else
    switch optimizer
      case 'GaussNewton'
        % Nothing to do here.
      case 'lsqnonlin'
        bounds = ones(length(params0),1)*[-Inf,Inf];
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),1)  = 0.0;
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),2)  = 1.0;
      case 'fmincon'
        if isoctave && ~user_has_octave_forge_package('optim', '1.6')
            error('Optimization algorithm ''fmincon'' requires the optim package, version 1.6 or higher')
        elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
            error('Optimization algorithm ''fmincon'' requires the Optimization Toolbox')
        end
        minalgo = 1;
        bounds = ones(length(params0),1)*[-Inf,Inf];
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),1)  = 0.0;
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),2)  = 1.0;
      case 'fminunc'
        if isoctave && ~user_has_octave_forge_package('optim')
            error('Optimization algorithm ''fminunc'' requires the optim package')
        elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
            error('Optimization algorithm ''fminunc'' requires the Optimization Toolbox')
        end
        minalgo = 3;
      case 'csminwel'
        minalgo = 4;
      case 'fminsearch'
        if isoctave && ~user_has_octave_forge_package('optim')
            error('Optimization algorithm ''fminsearch'' requires the optim package')
        elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
            error('Optimization algorithm ''fminsearch'' requires the Optimization Toolbox')
        end
        minalgo = 7;
      case 'simplex'
        minalgo = 8;
      case 'annealing'
        minalgo = 2;
        bounds = ones(length(params0),1)*[-Inf,Inf];
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),1)  = 0.0;
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),2)  = 1.0;
        parameter_names = fieldnames(params);
      case 'particleswarm'
        if isoctave
            error('Optimization ''particleswarm'' is not available under Octave')
        elseif ~user_has_matlab_license('global_optimization_toolbox')
            error('Optimization ''particleswarm'' requires the Global Optimization Toolbox')
        end
        minalgo = 12;
        bounds = ones(length(params0),1)*[-Inf,Inf];
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),1)  = 0.0;
        bounds(strcmp(fieldnames(params), M_.param_names(M_.pac.(pacmodl).ec.params)),1)  = 1.0;
        parameter_names = fieldnames(params);
      otherwise
        msg = sprintf('%s is not an implemented optimization routine.\n', optimizer);
        msg = sprintf('%sAvailable algorithms are:\n', msg);
        msg = sprintf('%s - %s\n', msg, 'fmincon');
        msg = sprintf('%s - %s\n', msg, 'fminunc');
        msg = sprintf('%s - %s\n', msg, 'csminwel');
        msg = sprintf('%s - %s\n', msg, 'fminsearch');
        msg = sprintf('%s - %s\n', msg, 'simplex');
        msg = sprintf('%s - %s\n', msg, 'annealing');
        msg = sprintf('%s - %s\n', msg, 'lsqnonlin');
        msg = sprintf('%s - %s\n', msg, 'GaussNewton');
        error(msg)
    end
end

% Set options if provided as input arguments to nls routine.
oldopt = options_.optim_opt;
[noprint, opt] = opt4nls(varargin);
options_.optim_opt = opt;

%
% Check that we are able to evaluate the Sum of Squared Residuals on the initial guess
%

ssr0 = ssrfun(params0);

if isnan(ssr0) || isinf(ssr0) || ~isreal(ssr0)
    error('Cannot evaluate the Sum of Squared Residuals on the initial guess.')
end

if is_gauss_newton
    [params1, SSR] = gauss_newton(resfun, params0);
elseif is_lsqnonlin
    if ismember('levenberg-marquardt', varargin)
        % Levenberg Marquardt does not handle boundary constraints.
        [params1, SSR] = lsqnonlin(resfun, params0, [], [], optimset(varargin{:}));
    else
        [params1, SSR] = lsqnonlin(resfun, params0, bounds(:,1), bounds(:,2), optimset(varargin{:}));
    end
else
    % Estimate the parameters by minimizing the sum of squared residuals.
    [params1, SSR] = dynare_minimize_objective(ssrfun, params0, ...
                                               minalgo, ...
                                               options_, ...
                                               bounds, ...
                                               parameter_names, ...
                                               [], ...
                                               []);
end

% Revert local modifications to the options.
options_.optim_opt = oldopt;

% Compute an estimator of the covariance matrix (see White and
% Domovitz [Econometrica, 1984], theorem 3.2).
[r, J] = jacobian(resfun, params1, 1e-6);
T = length(r);
A = 2.0*(J'*J)/T;
J = bsxfun(@times, J, r);
B = J'*J;
l = round(T^.25);
for tau=1:l
    B = B + (1-tau/(l+1))*(J(tau+1:end,:)'*J(1:end-tau,:)+J(1:end-tau,:)'*J(tau+1:end,:));
end
B = (4.0/T)*B;
C = inv(A)*B*inv(A); % C is the asymptotic covariance of sqrt(T) times the vector of estimated parameters.
C = C/T;

% Save results
lhs = eval(strrep(lhs, 'data', 'data(range(1)-1:range(end)).data'));
oo_.pac.(pacmodl).lhs = dseries(lhs, range(1), 'lhs');
oo_.pac.(pacmodl).fit = dseries(lhs-r, range(1), 'fit');
oo_.pac.(pacmodl).residual = dseries(r, range(1), 'residual');
oo_.pac.(pacmodl).ssr = SSR;
oo_.pac.(pacmodl).s2 = SSR/T;
oo_.pac.(pacmodl).R2 = 1-var(r)/var(lhs);
oo_.pac.(pacmodl).parnames = fieldnames(params);
oo_.pac.(pacmodl).estimator = params1;
oo_.pac.(pacmodl).covariance = C;
oo_.pac.(pacmodl).student = params1./(sqrt(diag(C)));

% Also save estimated parameters in M_
M_.params(ipnames_) = params1;
M_ = pac.update.parameters(pacmodl, M_, oo_, false);

if ~noprint
    title = ['NLS Estimation of equation ''' eqname ''''];
    preamble = {['Dependent Variable: ' dLHS], ...
                sprintf('Observations: %d from %s to %s\n', (range(end)-range(1))+1, range(1).char, range(end).char)};

    afterward = {sprintf('R^2: %f', oo_.pac.(pacmodl).R2), ...
                 sprintf('s^2: %f', oo_.pac.(pacmodl).s2)}; ...

    dyn_table(title, preamble, afterward, oo_.pac.(pacmodl).parnames, ...
              {'Estimates','t-statistic','Std. Error'}, 4, ...
              [oo_.pac.(pacmodl).estimator oo_.pac.(pacmodl).student sqrt(diag(C))]);
end
