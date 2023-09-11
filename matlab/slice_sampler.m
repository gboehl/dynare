function [theta, fxsim, neval] = slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)
% function [theta, fxsim, neval] = slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)
% ----------------------------------------------------------
% UNIVARIATE SLICE SAMPLER - stepping out (Neal, 2003)
% W: optimal value in the range (3,10)*std(x)
%    - see C.Planas and A.Rossi (2014)
% objective_function(theta,varargin): -log of any unnormalized pdf
% with varargin (optional) a vector of auxiliaty parameters
% to be passed to f( ).
% ----------------------------------------------------------
%
% INPUTS
%   objective_function:       objective function (expressed as minus the log of a density)
%   theta:                    last value of theta
%   thetaprior:               bounds of the theta space
%   sampler_options:          posterior sampler options
%   varargin:                 optional input arguments to objective function
%
% OUTPUTS
%   theta:       new theta sample
%   fxsim:       value of the objective function for the new sample
%   neval:       number of function evaluations
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2015-2023 Dynare Team
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

if sampler_options.rotated %&& ~isempty(sampler_options.V1),
    [theta, fxsim, neval] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin{:});
    if isempty(sampler_options.mode) % jumping
        return
    else
        nevalR=sum(neval);
    end
end

theta=theta(:);
npar = length(theta);
W1 = sampler_options.W1;
neval = zeros(npar,1);

% % % fname0=fname;
fname = [ int2str(sampler_options.curr_block)];

Prior = dprior(varargin{6},varargin{3}.prior_trunc);

it=0;
while it<npar
    it=it+1;
    neval(it) = 0;
    W = W1(it);
    xold  = theta(it);
    % XLB   = thetaprior(3);
    % XUB   = thetaprior(4);
    XLB   = thetaprior(it,1);
    XUB   = thetaprior(it,2);


    % -------------------------------------------------------
    % 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
    %    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
    % -------------------------------------------------------
    fxold = -feval(objective_function,theta,varargin{:});
    if ~isfinite(fxold)
        disp(['slice_sampler:: Iteration ' int2str(it) ' started with bad parameter set (fval is inf or nan)'])
        icount=0;
        while ~isfinite(fxold) && icount<1000
            icount=icount+1;
            theta = Prior.draw();
            if all(theta >= thetaprior(:,1)) && all(theta <= thetaprior(:,2))
                fxold = -feval(objective_function,theta,varargin{:});
            end
        end
        % restart from 1
        it = 1;
        neval(it) = 0;
        W = W1(it);
        xold  = theta(it);
        XLB   = thetaprior(it,1);
        XUB   = thetaprior(it,2);
    end
    neval(it) = neval(it) + 1;
    Z = fxold + log(rand(1,1));
    % -------------------------------------------------------------
    % 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE
    %    STEPPING-OUT PROCEDURE
    % -------------------------------------------------------------
    u = rand(1,1);
    L = max(XLB,xold-W*u);
    R = min(XUB,L+W);
    mytxt{it,1} = '';
    while(L > XLB)
        xsim = L;
        theta(it) = xsim;
        fxl = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (fxl <= Z)
            break
        end
        L = max(XLB,L-W);
        if neval(it)>30
            L=XLB;
            xsim = L;
            theta(it) = xsim;
            fxl = -feval(objective_function,theta,varargin{:});
            icount = 0;
            while (isinf(fxl) || isnan(fxl)) && icount<300
                icount = icount+1;
                L=L+sqrt(eps);
                xsim = L;
                theta(it) = xsim;
                fxl = -feval(objective_function,theta,varargin{:});
            end
            mytxt{it,1} = sprintf('Getting L for [%s] is taking too long.', varargin{6}.name{it});
            save(['slice_iter_info_' fname],'mytxt','neval','it','theta','fxl')
            %keyboard;
        end
    end
    neval1 = neval(it);
    mytxt{it,2} = '';
    while(R < XUB)
        xsim = R;
        theta(it) = xsim;
        fxr = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (fxr <= Z)
            break
        end
        R = min(XUB,R+W);
        if neval(it)>(neval1+30)
            R=XUB;
            xsim = R;
            theta(it) = xsim;
            fxr = -feval(objective_function,theta,varargin{:});
            icount = 0;
            while (isinf(fxr) || isnan(fxr)) && icount<300
                icount = icount+1;
                R=R-sqrt(eps);
                xsim = R;
                theta(it) = xsim;
                fxr = -feval(objective_function,theta,varargin{:});
            end
            mytxt{it,2} = sprintf('Getting R for [%s] is taking too long.', varargin{6}.name{it});
            save(['slice_iter_info_' fname],'mytxt','neval','it','theta','fxr')
            %keyboard;
        end
    end
    % ------------------------------------------------------
    % 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
    % ------------------------------------------------------
    fxsim = Z-1;
    mytxt{it,3} = '';
    while (fxsim < Z)
        u = rand(1,1);
        xsim = L + u*(R - L);
        theta(it) = xsim;
        fxsim = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (xsim > xold)
            R = xsim;
        else
            L = xsim;
        end
        if (R-L)<1.e-6 %neval(it)>(30+neval2)
            fprintf('The sampling for parameter [%s] is taking too long as the sampling set is too tight. Check the prior.\n', varargin{6}.name{it})
            mytxt{it,3} = sprintf('Sampling [%s] is taking too long.', varargin{6}.name{it});
            save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it')
            break
        end
    end
    save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxsim')

    if isinf(fxsim) || isnan(fxsim)
        theta(it) = xold;
        fxsim = fxold;
        disp('SLICE: posterior density is infinite. Reset values at initial ones.')
    end
end

if sampler_options.rotated && ~isempty(sampler_options.mode) % jumping
    neval=sum(neval)+nevalR;
end
