function Scale = calibrate_mh_scale_parameter(ObjectiveFunction, CovarianceMatrix, Parameters, MhBounds, options, varargin)

% Tune the MH scale parameter so that the overall acceptance ratio is close to AcceptanceTarget.
%
% INPUTS
% - ObjectiveFunction             [fhandle]      Function (posterior kernel).
% - CovarianceMatrix              [double]       n*n matrix, covariance matrix of the jumping distribution.
% - Parameters                    [double]       n*1 vector, parameter values.
% - MhBounds                      [double]       n*2 matrix, bounds on the possible values for the parameters.
% - options                       [structure]    content of options_.tune_mh_jscale.
% - varargin                      [cell]         Additional arguments to be passed to ObjectiveFunction.
%
% OUTPUTS
% - Scale                         [double]       scalar, optimal scale parameter for teh jumping distribution.

% Fire up the wait bar
hh = dyn_waitbar(0,'Tuning of the scale parameter...');
set(hh,'Name','Tuning of the scale parameter.');

% Intilialize various counters.
j = 1; jj  = 1; isux = 0; jsux = 0; i = 0;

% Evaluate the objective function.
logpo0 = - feval(ObjectiveFunction, Parameters, varargin{:});
logpo1 = logpo0;

% Get the dimension of the problem.
n = length(Parameters);

% Initialize the correction on the scale factor.
correction = 1.0;

% Set the initial value of the scale parameter
Scale = options.guess;

% Transposition of some arrays.
MhBounds = MhBounds';
Parameters = Parameters';

% Compute the Cholesky of the covariance matrix, return an error if the
% matrix is not positive definite.
try
    dd = chol(CovarianceMatrix);
catch
    error('The covariance matrix has to be a symetric positive definite matrix!')
end

% Set parameters related to the proposal distribution
if options.rwmh.proposal_distribution=='rand_multivariate_normal'
    nu = n;
elseif options.rwmh.proposal_distribution=='rand_multivariate_student'
    nu = options.rwmh.student_degrees_of_freedom;
end

% Random Walk Metropolis Hastings iterations...
while j<=options.maxiter
    % Obtain a proposal (jump)
    proposal = feval(options.rwmh.proposal_distribution, Parameters, Scale*dd, nu);
    % If out of boundaries set the posterior kernel equal to minus infinity
    % so that the proposal will be rejected with probability one.
    if all(proposal > MhBounds(1,:)) && all(proposal < MhBounds(2,:))
        logpo0 = -feval(ObjectiveFunction, proposal(:), varargin{:});
    else
        logpo0 = -inf;
    end
    % Move if the proposal is enough likely...
    if logpo0>-inf && log(rand)<logpo0-logpo1
        Parameters = proposal;
        logpo1 = logpo0;
        isux = isux + 1;
        jsux = jsux + 1;
    end% ... otherwise I don't move.
    prtfrc = j/options.maxiter;
    % Update the waitbar
    if ~mod(j, 10)
        dyn_waitbar(prtfrc, hh, sprintf('Acceptance ratio [during last 500]: %f [%f]', isux/j, jsux/jj));
    end
    % Adjust the value of the scale parameter.
    if ~mod(j, options.stepsize)
        r1 = jsux/jj;    % Local acceptance ratio
        r2 = isux/j;     % Overall acceptance ratio
        % Set correction for the scale factor
        c1 = r1/options.target;
        if abs(c1-1)>.05
            correction = correction^options.rho*c1^(1-options.rho);
        else
            correction = c1;
        end
        % Apply correction
        if c1>0
            Scale = Scale*correction;
        else
            Scale = Scale/10;
        end
        % Update some counters.
        jsux = 0; jj = 0;
        if abs(r2-options.target)<options.c2 && abs(r1-options.target)<options.c1
            i = i+1;
        else
            i = 0;
        end
        % Test convergence.
        if i>options.c3
            break
        end
    end
    j = j+1;
    jj = jj + 1;
end

dyn_waitbar_close(hh);