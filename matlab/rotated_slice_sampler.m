function [theta, fxsim, neval] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)
% ----------------------------------------------------------
% ROTATED SLICE SAMPLER - with stepping out (Neal, 2003)
% extension of the orthogonal univarite sampler (slice_sampler.m)
% copyright M. Ratto (European Commission)
%
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

theta=theta(:);
npar = length(theta);
neval = zeros(npar,1);
W1=[];
if isfield(sampler_options,'WR')
    W1 = sampler_options.WR;
end
if ~isempty(sampler_options.mode)
    mm = sampler_options.mode;
    n = length(mm);
    distance=NaN(n,1);
    for j=1:n
        distance(j)=sqrt(sum((theta-mm(j).m).^2));
    end
    [~, im] = min(distance);

    r=im;
    V1 = mm(r).m;
    jj=0;
    for j=1:n
        if j~=r
            jj=jj+1;
            tmp=mm(j).m-mm(r).m;
            %tmp=mm(j).m-theta;
            V1(:,jj)=tmp/norm(tmp);
        end
    end
    resul=randperm(n-1,n-1);
    V1 = V1(:,resul);
else
    V1 = sampler_options.V1;
end
npar=size(V1,2);

for it=1:npar
    theta0 = theta;
    neval(it) = 0;
    xold  = 0;
    tb=sort([(thetaprior(:,1)-theta)./V1(:,it) (thetaprior(:,2)-theta)./V1(:,it)],2);
    XLB=max(tb(:,1));
    XUB=min(tb(:,2));
    if isempty(W1)
        W = (XUB-XLB); %*0.8;
    else
        W = W1(it);
    end

    % -------------------------------------------------------
    % 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
    %    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
    % -------------------------------------------------------

    fxold = -feval(objective_function,theta,varargin{:});
    %I have to be sure that the rotation is for L,R or for Fxold, theta(it)
    neval(it) = neval(it) + 1;
    Z = fxold + log(rand(1,1));
    % -------------------------------------------------------------
    % 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE
    %    STEPPING-OUT PROCEDURE
    % -------------------------------------------------------------
    u = rand(1,1);
    L = max(XLB,xold-W*u);
    R = min(XUB,L+W);

    %[L R]=slice_rotation(L, R, alpha);
    while(L > XLB)
        xsim = L;
        theta = theta0+xsim*V1(:,it);
        fxl = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (fxl <= Z)
            break
        end
        L = max(XLB,L-W);
    end
    while(R < XUB)
        xsim = R;
        theta = theta0+xsim*V1(:,it);
        fxr = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (fxr <= Z)
            break
        end
        R = min(XUB,R+W);
    end
    % ------------------------------------------------------
    % 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
    % ------------------------------------------------------
    fxsim = Z-1;
    while (fxsim < Z)
        u = rand(1,1);
        xsim = L + u*(R - L);
        theta = theta0+xsim*V1(:,it);
        fxsim = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (xsim > xold)
            R = xsim;
        else
            L = xsim;
        end
    end
end
end
