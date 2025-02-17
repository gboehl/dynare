function [X,FVAL,EXITFLAG,OUTPUT] = simpsa(FUN,X0,LB,UB,OPTIONS,varargin)

% Finds a minimum of a function of several variables using an algorithm
% that is based on the combination of the non-linear smplex and the simulated
% annealing algorithm (the SIMPSA algorithm, Cardoso et al., 1996).
% In this paper, the algorithm is shown to be adequate for the global optimi-
% zation of an example set of unconstrained and constrained NLP functions.
%
%   SIMPSA attempts to solve problems of the form:
%       min F(X) subject to: LB <= X <= UB
%        X
%
% Algorithm partly is based on paper of Cardoso et al, 1996.
%
%   X=SIMPSA(FUN,X0) start at X0 and finds a minimum X to the function FUN.
%   FUN accepts input X and returns a scalar function value F evaluated at X.
%   X0 may be a scalar, vector, or matrix.
%
%   X=SIMPSA(FUN,X0,LB,UB) defines a set of lower and upper bounds on the
%   design variables, X, so that a solution is found in the range
%   LB <= X <= UB. Use empty matrices for LB and UB if no bounds exist.
%   Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if X(i) is
%   unbounded above.
%
%   X=SIMPSA(FUN,X0,LB,UB,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the SIMPSASET function. See SIMPSASET for details.
%   Used options are TEMP_START, TEMP_END, COOL_RATE, INITIAL_ACCEPTANCE_RATIO,
%   MIN_COOLING_FACTOR, MAX_ITER_TEMP_FIRST, MAX_ITER_TEMP_LAST, MAX_ITER_TEMP,
%   MAX_ITER_TOTAL, MAX_TIME, MAX_FUN_EVALS, TOLX, TOLFUN, DISPLAY and OUTPUT_FCN.
%   Use OPTIONS = [] as a place holder if no options are set.
%
%   X=SIMPSA(FUN,X0,LB,UB,OPTIONS,varargin) is used to supply a variable
%   number of input arguments to the objective function FUN.
%
%   [X,FVAL]=SIMPSA(FUN,X0,...) returns the value of the objective
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG]=SIMPSA(FUN,X0,...) returns an EXITFLAG that describes the
%   exit condition of SIMPSA. Possible values of EXITFLAG and the corresponding
%   exit conditions are:
%
%     1  Change in the objective function value less than the specified tolerance.
%     2  Change in X less than the specified tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Maximum time exceeded.
%
%   [X,FVAL,EXITFLAG,OUTPUT]=SIMPSA(FUN,X0,...) returns a structure OUTPUT with
%   the number of iterations taken in OUTPUT.nITERATIONS, the number of function
%   evaluations in OUTPUT.nFUN_EVALS, the temperature profile in OUTPUT.TEMPERATURE,
%   the simplexes that were evaluated in OUTPUT.SIMPLEX and the best one in
%   OUTPUT.SIMPLEX_BEST, the costs associated with each simplex in OUTPUT.COSTS and
%   from the best simplex at that iteration in OUTPUT.COST_BEST, the amount of time
%   needed in OUTPUT.TIME and the options used in OUTPUT.OPTIONS.
%
%   See also SIMPSASET, SIMPSAGET


% Copyright © 2005 Henning Schmidt, FCC, henning@fcc.chalmers.se
% Copyright © 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
% Copyright © 2013-2017 Dynare Team.
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

% handle variable input arguments

if nargin < 5
    OPTIONS = [];
    if nargin < 4
        UB = 1e5;
        if nargin < 3
            LB = -1e5;
        end
    end
end

% check input arguments

if ~ischar(FUN)
    error('''FUN'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(X0)
    error('''X0'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(LB)
    error('''LB'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(UB)
    error('''UB'' incorrectly specified in ''SIMPSA''');
end
if length(X0) ~= length(LB)
    error('''LB'' and ''X0'' have incompatible dimensions in ''SIMPSA''');
end
if length(X0) ~= length(UB)
    error('''UB'' and ''X0'' have incompatible dimensions in ''SIMPSA''');
end

% declaration of global variables

global NDIM nFUN_EVALS TEMP YBEST PBEST

% set EXITFLAG to default value

EXITFLAG = -2;

% determine number of variables to be optimized

NDIM = length(X0);

% set default options
DEFAULT_OPTIONS = simpsaset('TEMP_START',[],...  % starting temperature (if none provided, an optimal one will be estimated)
                            'TEMP_END',.1,...                    % end temperature
                            'COOL_RATE',10,...                  % small values (<1) means slow convergence,large values (>1) means fast convergence
                            'INITIAL_ACCEPTANCE_RATIO',0.95,... % when initial temperature is estimated, this will be the initial acceptance ratio in the first round
                            'MIN_COOLING_FACTOR',0.9,...        % minimum cooling factor (<1)
                            'MAX_ITER_TEMP_FIRST',50,...        % number of iterations in the preliminary temperature loop
                            'MAX_ITER_TEMP_LAST',2000,...         % number of iterations in the last temperature loop (pure simplex)
                            'MAX_ITER_TEMP',10,...              % number of iterations in the remaining temperature loops
                            'MAX_ITER_TOTAL',2500,...           % maximum number of iterations tout court
                            'MAX_TIME',2500,...                 % maximum duration of optimization
                            'MAX_FUN_EVALS',20000,...            % maximum number of function evaluations
                            'TOLX',1e-6,...                     % maximum difference between best and worst function evaluation in simplex
                            'TOLFUN',1e-6,...                   % maximum difference between the coordinates of the vertices
                            'DISPLAY','iter',...                % 'iter' or 'none' indicating whether user wants feedback
                            'OUTPUT_FCN',[]);                   % string with output function name

% update default options with supplied options

OPTIONS = simpsaset(DEFAULT_OPTIONS,OPTIONS);

% store options in OUTPUT
if nargout>3
    OUTPUT.OPTIONS = OPTIONS;
end

% initialize simplex
% ------------------

% create empty simplex matrix p (location of vertex i in row i)
P = zeros(NDIM+1,NDIM);
% create empty cost vector (cost of vertex i in row i)
Y = zeros(NDIM+1,1);
% set best vertex of initial simplex equal to initial parameter guess
PBEST = X0(:)';
% calculate cost with best vertex of initial simplex
YBEST = CALCULATE_COST(FUN,PBEST,LB,UB,varargin{:});

% initialize temperature loop
% ---------------------------

% set temperature loop number to one
TEMP_LOOP_NUMBER = 1;

% if no TEMP_START is supplied, the initial temperature is estimated in the first
% loop as described by Cardoso et al., 1996 (recommended)

% therefore, the temperature is set to YBEST*1e5 in the first loop
if isempty(OPTIONS.TEMP_START)
    TEMP = abs(YBEST)*1e5;
else
    TEMP = OPTIONS.TEMP_START;
end

% initialize OUTPUT structure
% ---------------------------
if nargout>3
    OUTPUT.TEMPERATURE = zeros(OPTIONS.MAX_ITER_TOTAL,1);
    OUTPUT.SIMPLEX = zeros(NDIM+1,NDIM,OPTIONS.MAX_ITER_TOTAL);
    OUTPUT.SIMPLEX_BEST = zeros(OPTIONS.MAX_ITER_TOTAL,NDIM);
    OUTPUT.COSTS = zeros(OPTIONS.MAX_ITER_TOTAL,NDIM+1);
    OUTPUT.COST_BEST = zeros(OPTIONS.MAX_ITER_TOTAL,1);
end
% initialize iteration data
% -------------------------

% start timer
tic
% set number of function evaluations to one
nFUN_EVALS = 1;
% set number of iterations to zero
nITERATIONS = 0;

% temperature loop: run SIMPSA till stopping criterion is met
% -----------------------------------------------------------

while 1

    % detect if termination criterium was met
    % ---------------------------------------

    % if a termination criterium was met, the value of EXITFLAG should have changed
    % from its default value of -2 to -1, 0, 1 or 2

    if EXITFLAG ~= -2
        break
    end

    % set MAXITERTEMP: maximum number of iterations at current temperature
    % --------------------------------------------------------------------

    if TEMP_LOOP_NUMBER == 1
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP_FIRST*NDIM;
        % The initial temperature is estimated (is requested) as described in
        % Cardoso et al. (1996). Therefore, we need to store the number of
        % successful and unsuccessful moves, as well as the increase in cost
        % for the unsuccessful moves.
        if isempty(OPTIONS.TEMP_START)
            [SUCCESSFUL_MOVES,UNSUCCESSFUL_MOVES,UNSUCCESSFUL_COSTS] = deal(0);
        end
    elseif TEMP < OPTIONS.TEMP_END
        TEMP = 0;
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP_LAST*NDIM;
    else
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP*NDIM;
    end

    % construct initial simplex
    % -------------------------

    % 1st vertex of initial simplex
    P(1,:) = PBEST;
    Y(1) = CALCULATE_COST(FUN,P(1,:),LB,UB,varargin{:});

    % if output function given then run output function to plot intermediate result
    if ~isempty(OPTIONS.OUTPUT_FCN)
        feval(OPTIONS.OUTPUT_FCN,transpose(P(1,:)),Y(1));
    end

    % remaining vertices of simplex
    for k = 1:NDIM
        % copy first vertex in new vertex
        P(k+1,:) = P(1,:);
        % alter new vertex
        P(k+1,k) = LB(k)+rand*(UB(k)-LB(k));
        % calculate value of objective function at new vertex
        Y(k+1) = CALCULATE_COST(FUN,P(k+1,:),LB,UB,varargin{:});
    end

    % store information on what step the algorithm just did
    ALGOSTEP = 'initial simplex';

    % add NDIM+1 to number of function evaluations
    nFUN_EVALS = nFUN_EVALS + NDIM;

    % note:
    %  dimensions of matrix P: (NDIM+1) x NDIM
    %  dimensions of vector Y: (NDIM+1) x 1

    % give user feedback if requested
    if strcmp(OPTIONS.DISPLAY,'iter')
        if nITERATIONS == 0
            disp(' Nr Iter  Nr Fun Eval    Min function       Best function        TEMP           Algorithm Step');
        else
            dprintf('%5.0f      %5.0f       %12.6g     %15.6g      %12.6g       %s',nITERATIONS,nFUN_EVALS,Y(1),YBEST,TEMP,'best point');
        end
    end

    % run full metropolis cycle at current temperature
    % ------------------------------------------------

    % initialize vector COSTS, needed to calculate new temperature using cooling
    % schedule as described by Cardoso et al. (1996)
    COSTS = zeros((NDIM+1)*MAXITERTEMP,1);

    % initialize ITERTEMP to zero

    ITERTEMP = 0;

    % start

    for ITERTEMP = 1:MAXITERTEMP

        % add one to number of iterations
        nITERATIONS = nITERATIONS + 1;

        % Press and Teukolsky (1991) add a positive logarithmic distributed variable,
        % proportional to the control temperature T to the function value associated with
        % every vertex of the simplex. Likewise,they subtract a similar random variable
        % from the function value at every new replacement point.
        % Thus, if the replacement point corresponds to a lower cost, this method always
        % accepts a true down hill step. If, on the other hand, the replacement point
        % corresponds to a higher cost, an uphill move may be accepted, depending on the
        % relative COSTS of the perturbed values.
        % (taken from Cardoso et al.,1996)

        % add random fluctuations to function values of current vertices
        YFLUCT = Y+TEMP*abs(log(rand(NDIM+1,1)));

        % reorder YFLUCT, Y and P so that the first row corresponds to the lowest YFLUCT value
        help = sortrows([YFLUCT,Y,P],1);
        YFLUCT = help(:,1);
        Y = help(:,2);
        P = help(:,3:end);

        if nargout>3
            % store temperature at current iteration
            OUTPUT.TEMPERATURE(nITERATIONS) = TEMP;

            % store information about simplex at the current iteration
            OUTPUT.SIMPLEX(:,:,nITERATIONS) = P;
            OUTPUT.SIMPLEX_BEST(nITERATIONS,:) = PBEST;

            % store cost function value of best vertex in current iteration
            OUTPUT.COSTS(nITERATIONS,:) = Y;
            OUTPUT.COST_BEST(nITERATIONS) = YBEST;
        end

        if strcmp(OPTIONS.DISPLAY,'iter')
            dprintf('%5.0f      %5.0f       %12.6g     %15.6g      %12.6g       %s',nITERATIONS,nFUN_EVALS,Y(1),YBEST,TEMP,ALGOSTEP);
        end

        % if output function given then run output function to plot intermediate result
        if ~isempty(OPTIONS.OUTPUT_FCN)
            feval(OPTIONS.OUTPUT_FCN,transpose(P(1,:)),Y(1));
        end

        % end the optimization if one of the stopping criteria is met
        %% 1. difference between best and worst function evaluation in simplex is smaller than TOLFUN
        %% 2. maximum difference between the coordinates of the vertices in simplex is less than TOLX
        %% 3. no convergence,but maximum number of iterations has been reached
        %% 4. no convergence,but maximum time has been reached

        if (abs(max(Y)-min(Y)) < OPTIONS.TOLFUN) && (TEMP_LOOP_NUMBER ~= 1)
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Change in the objective function value less than the specified tolerance (TOLFUN).')
            end
            EXITFLAG = 1;
            break
        end

        if (max(max(abs(P(2:NDIM+1,:)-P(1:NDIM,:)))) < OPTIONS.TOLX) && (TEMP_LOOP_NUMBER ~= 1)
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Change in X less than the specified tolerance (TOLX).')
            end
            EXITFLAG = 2;
            break
        end

        if (nITERATIONS >= OPTIONS.MAX_ITER_TOTAL*NDIM) || (nFUN_EVALS >= OPTIONS.MAX_FUN_EVALS*NDIM*(NDIM+1))
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Maximum number of function evaluations or iterations reached.');
            end
            EXITFLAG = 0;
            break
        end

        if toc/60 > OPTIONS.MAX_TIME
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Exceeded maximum time.');
            end
            EXITFLAG = -1;
            break
        end

        % begin a new iteration

        %% first extrapolate by a factor -1 through the face of the simplex
        %% across from the high point,i.e.,reflect the simplex from the high point
        [YFTRY,YTRY,PTRY] = AMOTRY(FUN,P,-1,LB,UB,varargin{:});

        %% check the result
        if YFTRY <= YFLUCT(1)
            %% gives a result better than the best point,so try an additional
            %% extrapolation by a factor 2
            [YFTRYEXP,YTRYEXP,PTRYEXP] = AMOTRY(FUN,P,-2,LB,UB,varargin{:});
            if YFTRYEXP < YFTRY
                P(end,:) = PTRYEXP;
                Y(end) = YTRYEXP;
                ALGOSTEP = 'reflection and expansion';
            else
                P(end,:) = PTRY;
                Y(end) = YTRY;
                ALGOSTEP = 'reflection';
            end
        elseif YFTRY >= YFLUCT(NDIM)
            %% the reflected point is worse than the second-highest, so look
            %% for an intermediate lower point, i.e., do a one-dimensional
            %% contraction
            [YFTRYCONTR,YTRYCONTR,PTRYCONTR] = AMOTRY(FUN,P,-0.5,LB,UB,varargin{:});
            if YFTRYCONTR < YFLUCT(end)
                P(end,:) = PTRYCONTR;
                Y(end) = YTRYCONTR;
                ALGOSTEP = 'one dimensional contraction';
            else
                %% can't seem to get rid of that high point, so better contract
                %% around the lowest (best) point
                X = ones(NDIM,NDIM)*diag(P(1,:));
                P(2:end,:) = 0.5*(P(2:end,:)+X);
                for k=2:NDIM
                    Y(k) = CALCULATE_COST(FUN,P(k,:),LB,UB,varargin{:});
                end
                ALGOSTEP = 'multiple contraction';
            end
        else
            %% if YTRY better than second-highest point, use this point
            P(end,:) = PTRY;
            Y(end) = YTRY;
            ALGOSTEP = 'reflection';
        end

        % the initial temperature is estimated in the first loop from
        % the number of successfull and unsuccesfull moves, and the average
        % increase in cost associated with the unsuccessful moves

        if TEMP_LOOP_NUMBER == 1 && isempty(OPTIONS.TEMP_START)
            if Y(1) > Y(end)
                SUCCESSFUL_MOVES = SUCCESSFUL_MOVES+1;
            elseif Y(1) <= Y(end)
                UNSUCCESSFUL_MOVES = UNSUCCESSFUL_MOVES+1;
                UNSUCCESSFUL_COSTS = UNSUCCESSFUL_COSTS+(Y(end)-Y(1));
            end
        end

    end

    % stop if previous for loop was broken due to some stop criterion
    if ITERTEMP < MAXITERTEMP
        break
    end

    % store cost function values in COSTS vector
    COSTS((ITERTEMP-1)*NDIM+1:ITERTEMP*NDIM+1) = Y;

    % calculated initial temperature or recalculate temperature
    % using cooling schedule as proposed by Cardoso et al. (1996)
    % -----------------------------------------------------------

    if TEMP_LOOP_NUMBER == 1 && isempty(OPTIONS.TEMP_START)
        TEMP = -(UNSUCCESSFUL_COSTS/(SUCCESSFUL_MOVES+UNSUCCESSFUL_MOVES))/log(((SUCCESSFUL_MOVES+UNSUCCESSFUL_MOVES)*OPTIONS.INITIAL_ACCEPTANCE_RATIO-SUCCESSFUL_MOVES)/UNSUCCESSFUL_MOVES);
    elseif TEMP_LOOP_NUMBER ~= 0
        STDEV_Y = std(COSTS);
        COOLING_FACTOR = 1/(1+TEMP*log(1+OPTIONS.COOL_RATE)/(3*STDEV_Y));
        TEMP = TEMP*min(OPTIONS.MIN_COOLING_FACTOR,COOLING_FACTOR);
    end

    % add one to temperature loop number
    TEMP_LOOP_NUMBER = TEMP_LOOP_NUMBER+1;

end

% return solution
X = transpose(PBEST);
FVAL = YBEST;

if nargout>3
    % store number of function evaluations
    OUTPUT.nFUN_EVALS = nFUN_EVALS;

    % store number of iterations
    OUTPUT.nITERATIONS = nITERATIONS;

    % trim OUTPUT data structure
    OUTPUT.TEMPERATURE(nITERATIONS+1:end) = [];
    OUTPUT.SIMPLEX(:,:,nITERATIONS+1:end) = [];
    OUTPUT.SIMPLEX_BEST(nITERATIONS+1:end,:) = [];
    OUTPUT.COSTS(nITERATIONS+1:end,:) = [];
    OUTPUT.COST_BEST(nITERATIONS+1:end) = [];

    % store the amount of time needed in OUTPUT data structure
    OUTPUT.TIME = toc;
end

return

% ==============================================================================

% AMOTRY FUNCTION
% ---------------

function [YFTRY,YTRY,PTRY] = AMOTRY(FUN,P,fac,LB,UB,varargin)
% Extrapolates by a factor fac through the face of the simplex across from
% the high point, tries it, and replaces the high point if the new point is
% better.

global NDIM TEMP

% calculate coordinates of new vertex
psum = sum(P(1:NDIM,:))/NDIM;
PTRY = psum*(1-fac)+P(end,:)*fac;

% evaluate the function at the trial point.
YTRY = CALCULATE_COST(FUN,PTRY,LB,UB,varargin{:});
% substract random fluctuations to function values of current vertices
YFTRY = YTRY-TEMP*abs(log(rand(1)));

return

% ==============================================================================

% COST FUNCTION EVALUATION
% ------------------------

function [YTRY] = CALCULATE_COST(FUN,PTRY,LB,UB,varargin)

global YBEST PBEST NDIM nFUN_EVALS

for i = 1:NDIM
    % check lower bounds
    if PTRY(i) < LB(i)
        YTRY = 1e12+(LB(i)-PTRY(i))*1e6;
        return
    end
    % check upper bounds
    if PTRY(i) > UB(i)
        YTRY = 1e12+(PTRY(i)-UB(i))*1e6;
        return
    end
end

% calculate cost associated with PTRY
YTRY = feval(FUN,PTRY(:),varargin{:});

% add one to number of function evaluations
nFUN_EVALS = nFUN_EVALS + 1;

% save the best point ever
if YTRY < YBEST
    YBEST = YTRY;
    PBEST = PTRY;
end

return
