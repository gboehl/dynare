function [x,f,exitflag,n_f_evals,n_grad_evals,n_constraint_evals,n_constraint_gradient_evals]=solvopt(x,fun,grad,func,gradc,optim,varargin)
% [x,f,options]=solvopt(x,fun,grad,func,gradc,options,varargin)
%
% The function SOLVOPT, developed by Alexei Kuntsevich and Franz Kappe,
% performs a modified version of Shor's r-algorithm in
% order to find a local minimum resp. maximum of a nonlinear function
% defined on the n-dimensional Euclidean space or % a solution of a nonlinear
% constrained problem:
% min { f(x): g(x) (<)= 0, g(x) in R(m), x in R(n) }
%
% Inputs:
%   x       n-vector (row or column) of the coordinates of the starting
%           point,
% fun       name of an M-file (M-function) which computes the value
%           of the objective function <fun> at a point x,
%           synopsis: f=fun(x)
% grad      indicator whether objective function provides the gradient
%           vector of the function <fun> at a point x,
% func      name of an M-file (M-function) which computes the MAXIMAL
%           RESIDUAL(!) for a set of constraints at a point x,
%           synopsis: fc=func(x)
% gradc     name of an M-file (M-function) which computes the gradient
%           vector for the maximal residual constraint at a point x,
%           synopsis: gc=gradc(x)
% optim     Options structure with fields:
%    optim.minimizer_indicator= H, where sign(H)=-1 resp. sign(H)=+1 means minimize
%        resp. maximize <fun> (valid only for unconstrained problem)
%        and H itself is a factor for the initial trial step size
%        (optim.minimizer_indicator=-1 by default),
%    optim.TolX= relative error for the argument in terms of the
%        infinity-norm (1.e-4 by default),
%    optim.TolFun= relative error for the function value (1.e-6 by default),
%    optim.MaxIter= limit for the number of iterations (15000 by default),
%    optim.verbosity= control of the display of intermediate results and error
%        resp. warning messages (default value is 0, i.e., no intermediate
%        output but error and warning messages, see more in the manual),
%    optim.TolXConstraint= admissible maximal residual for a set of constraints
%        (optim.TolXConstraint=1e-8 by default, see more in the manual),
%   *optim.SpaceDilation= the coefficient of space dilation (2.5 by default),
%   *optim.LBGradientStep= lower bound for the stepsize used for the difference
%        approximation of gradients (1e-12 by default, see more in the manual).
%  (* ... changes should be done with care)
%
% Outputs:
% x                     optimal parameter vector (row resp. column),
% f                     optimum function value
% exitflag:             the number of iterations, if positive,
%                       or an abnormal stop code, if negative (see more in the manual),
% n_f_evals:            number of objective evaluations
% n_grad_evals:         number of gradient evaluations,
% n_constraint_evals:   number of constraint function evaluations,
% n_constraint_gradient_evals   number of constraint gradient evaluations.
%
%
% Algorithm: Kuntsevich, A.V., Kappel, F., SolvOpt - The solver for local nonlinear optimization problems
% (version 1.1, Matlab, C, FORTRAN). University of Graz, Graz, 1997.
%
%
% Copyright © 1997-2008, Alexei Kuntsevich and Franz Kappel
% Copyright © 2008-2015 Giovanni Lombardo
% Copyright © 2015-2017 Dynare Team
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



% strings: ----{
errmes='SolvOpt error:';
wrnmes='SolvOpt warning:';
error1='No function name and/or starting point passed to the function.';
error2='Argument X has to be a row or column vector of dimension > 1.';
error30='<fun> returns an empty string.';
error31='Function value does not exist (NaN is returned).';
error32='Function equals infinity at the point.';
error40='<grad> returns an improper matrix. Check the dimension.';
error41='Gradient does not exist (NaN is returned by <grad>).';
error42='Gradient equals infinity at the starting point.';
error43='Gradient equals zero at the starting point.';
error50='<func> returns an empty string.';
error51='<func> returns NaN at the point.';
error52='<func> returns infinite value at the point.';
error60='<gradc> returns an improper vector. Check the dimension';
error61='<gradc> returns NaN at the point.';
error62='<gradc> returns infinite vector at the point.';
error63='<gradc> returns zero vector at an infeasible point.';
error5='Function is unbounded.';
error6='Choose another starting point.';
warn1= 'Gradient is zero at the point, but stopping criteria are not fulfilled.';
warn20='Normal re-setting of a transformation matrix.' ;
warn21='Re-setting due to the use of a new penalty coefficient.' ;
warn4= 'Iterations limit exceeded.';
warn31='The function is flat in certain directions.';
warn32='Trying to recover by shifting insensitive variables.';
warn09='Re-run from recorded point.';
warn08='Ravine with a flat bottom is detected.';
termwarn0='SolvOpt: Normal termination.';
termwarn1='SolvOpt: Termination warning:';
appwarn='The above warning may be reasoned by inaccurate gradient approximation';
endwarn=[...
    'Premature stop is possible. Try to re-run the routine from the obtained point.               ';...
    'Result may not provide the optimum. The function apparently has many extremum points.        ';...
    'Result may be inaccurate in the coordinates. The function is flat at the optimum.            ';...
    'Result may be inaccurate in a function value. The function is extremely steep at the optimum.'];
% ----}

% ARGUMENTS PASSED ----{
if nargin<2           % Function and/or starting point are not specified
    exitflag=-1;
    disp(errmes);
    disp(error1);
    return
end
if nargin<3
    app=1;             % No user-supplied gradients
elseif isempty(grad) || grad==0
    app=1;
else
    app=0;                     % Exact gradients are supplied
end

% OPTIONS ----{
doptions.minimizer_indicator=1;
doptions.TolX=1e-6; %accuracy of argument
doptions.TolFun=1e-6; %accuracy of function (see Solvopt p.29)
doptions.MaxIter=15000;
doptions.verbosity=1;
doptions.TolXConstraint=1.e-8;
doptions.SpaceDilation=2.5;
doptions.LBGradientStep=1.e-11;

if nargin<4
    optim=doptions;
elseif isempty(optim)
    optim=doptions;
end
% Check the values:
optim.TolX=max(optim.TolX,1.e-12);
optim.TolFun=max(optim.TolFun,1.e-12);
optim.TolX=max(optim.LBGradientStep*1.e2,optim.TolX);
optim.TolX=min(optim.TolX,1);
optim.TolFun=min(optim.TolFun,1);
optim.TolXConstraint=max(optim.TolXConstraint,1e-12);
optim.SpaceDilation=max([optim.SpaceDilation,1.5]);
optim.LBGradientStep=max(optim.LBGradientStep,1e-11);
% ----}

if isempty(func)
    constr=0;
else
    constr=1;                  % Constrained problem
    if isempty(gradc)
        appconstr=1;
    else
        appconstr=0;            % Exact gradients of constraints are supplied
    end
end
% ----}

% STARTING POINT ----{
if max(size(x))<=1
    disp(errmes);
    disp(error2);
    exitflag=-2;
    return
elseif size(x,2)==1
    n=size(x,1);
    x=x';
    trx=1;
elseif size(x,1)==1
    n=size(x,2);
    trx=0;
else
    disp(errmes);
    disp(error2);
    exitflag=-2;
    return
end
% ----}

% WORKING CONSTANTS AND COUNTERS ----{

n_f_evals=0; n_grad_evals=0;      % function and gradient calculations
if constr
    n_constraint_evals=0;
    n_constraint_gradient_evals=0;      % same for constraints
end
epsnorm=1.e-15;
epsnorm2=1.e-30;    % epsilon & epsilon^2

if constr, h1=-1;                  % NLP: restricted to minimization
    cnteps=optim.TolXConstraint;                % Max. admissible residual
else
    h1=sign(optim.minimizer_indicator);         % Minimize resp. maximize a function
end

k=0;                               % Iteration counter

wdef=1/optim.SpaceDilation-1;               % Default space transf. coeff.

%Gamma control ---{
ajb=1+.1/n^2;                    % Base I
ajp=20;
ajpp=ajp;                        % Start value for the power
ajs=1.15;                        % Base II
knorms=0; gnorms=zeros(1,10);    % Gradient norms stored
                                 %---}

%Display control ---{
if optim.verbosity<=0, dispdata=0;
    if optim.verbosity==-1
        dispwarn=0;
    else
        dispwarn=1;
    end
else
    dispdata=round(optim.verbosity);
    dispwarn=1;
end
ld=dispdata;
%---}

%Stepsize control ---{
dq=5.1;                          % Step divider (at f_{i+1}>gamma*f_{i})
du20=2;du10=1.5;du03=1.05;       % Step multipliers (at certain steps made)
kstore=3;nsteps=zeros(1,kstore); % Steps made at the last 'kstore' iterations
if app
    des=6.3;                 % Desired number of steps per 1-D search
else
    des=3.3;
end
mxtc=3;                          % Number of trial cycles (steep wall detect)
                                 %---}
termx=0; limxterm=50;              % Counter and limit for x-criterion

ddx   =max(1e-11,optim.LBGradientStep);      % stepsize for gradient approximation

low_bound=-1+1e-4;                 % Lower bound cosine used to detect a ravine

ZeroGrad=n*1.e-16;                 % Lower bound for a gradient norm

nzero=0;                           % Zero-gradient events counter
                                   % Lower bound for values of variables taking into account
lowxbound=max([optim.TolX,1e-3]);
% Lower bound for function values to be considered as making difference
lowfbound=optim.TolFun^2;
krerun=0;                          % Re-run events counter
detfr=optim.TolFun*100;              % relative error for f/f_{record}
detxr=optim.TolX*10;               % relative error for norm(x)/norm(x_{record})

warnno=0;                          % the number of warn.mess. to end with

kflat=0;                           % counter for points of flatness
stepvanish=0;                      % counter for vanished steps
stopf=0;
% ----}  End of setting constants
% ----}  End of the preamble

% COMPUTE THE FUNCTION  ( FIRST TIME ) ----{
if trx
    if app
        f=feval(fun,x',varargin{:});
    else
        [f,g]=feval(fun,x',varargin{:});
        n_grad_evals=n_grad_evals+1;
    end    
else
    if app
        f=feval(fun,x,varargin{:});
    else
        [f,g]=feval(fun,x,varargin{:});
        n_grad_evals=n_grad_evals+1;
    end
end
n_f_evals=n_f_evals+1;
if isempty(f)
    if dispwarn
        disp(errmes)
        disp(error30)
    end
    exitflag=-3;
    if trx
        x=x';
    end
    return
elseif isnan(f)
    if dispwarn
        disp(errmes)
        disp(error31)
        disp(error6)
    end
    exitflag=-3;
    if trx
        x=x';
    end
    return
elseif abs(f)==Inf
    if dispwarn
        disp(errmes)
        disp(error32)
        disp(error6)
    end
    exitflag=-3;
    if trx
        x=x';
    end
    return
end
xrec=x; frec=f;     % record point and function value
                    % Constrained problem
if constr,  fp=f; kless=0;
    if trx
        fc=feval(func,x');
    else
        fc=feval(func,x);
    end
    if isempty(fc)
        if dispwarn
            disp(errmes)
            disp(error50)
        end
        exitflag=-5;
        if trx
            x=x';
        end
        return
    elseif isnan(fc)
        if dispwarn
            disp(errmes)
            disp(error51)
            disp(error6)
        end
        exitflag=-5;
        if trx
            x=x';
        end
        return
    elseif abs(fc)==Inf
        if dispwarn
            disp(errmes)
            disp(error52)
            disp(error6)
        end
        exitflag=-5;
        if trx
            x=x';
        end
        return
    end
    n_constraint_evals=n_constraint_evals+1;
    PenCoef=1;                              % first rough approximation
    if fc<=cnteps
        FP=1; fc=0;             % feasible point
    else
        FP=0;                   % infeasible point
    end
    f=f+PenCoef*fc;
end
% ----}
% COMPUTE THE GRADIENT ( FIRST TIME ) ----{
if app
    deltax=h1*ddx*ones(size(x));
    if constr
        if trx
            g=apprgrdn(x',fp,fun,deltax',1,varargin{:});
        else
            g=apprgrdn(x ,fp,fun,deltax,1,varargin{:});
        end
    else
        if trx
            g=apprgrdn(x',f,fun,deltax',1,varargin{:});
        else
            g=apprgrdn(x ,f,fun,deltax,1,varargin{:});
        end
    end
    n_f_evals=n_f_evals+n;
else
    %done above
end
if size(g,2)==1, g=g'; end
ng=norm(g);
if size(g,2)~=n
    if dispwarn
        disp(errmes)
        disp(error40)
    end
    exitflag=-4;
    if trx
        x=x';
    end
    return
elseif isnan(ng)
    if dispwarn
        disp(errmes)
        disp(error41)
        disp(error6)
    end
    exitflag=-4;
    if trx
        x=x';
    end
    return
elseif ng==Inf
    if dispwarn
        disp(errmes)
        disp(error42)
        disp(error6)
    end
    exitflag=-4;
    if trx
        x=x';
    end
    return
elseif ng<ZeroGrad
    if dispwarn
        disp(errmes)
        disp(error43)
        disp(error6)
    end
    exitflag=-4;
    if trx
        x=x';
    end
    return
end
if constr
    if ~FP
        if appconstr
            deltax=sign(x); idx=find(deltax==0);
            deltax(idx)=ones(size(idx));
            deltax=ddx*deltax;
            if trx
                gc=apprgrdn(x',fc,func,deltax',0);
            else
                gc=apprgrdn(x ,fc,func,deltax ,0);
            end
            n_constraint_evals=n_constraint_evals+n;
        else
            if trx
                gc=feval(gradc,x');
            else
                gc=feval(gradc,x);
            end
            n_constraint_gradient_evals=n_constraint_gradient_evals+1;
        end
        if size(gc,2)==1
            gc=gc';
        end
        ngc=norm(gc);
        if size(gc,2)~=n
            if dispwarn
                disp(errmes)
                disp(error60)
            end
            exitflag=-6;
            if trx
                x=x';
            end
            return
        elseif isnan(ngc)
            if dispwarn
                disp(errmes)
                disp(error61)
                disp(error6)
            end
            exitflag=-6;
            if trx
                x=x';
            end
            return
        elseif ngc==Inf
            if dispwarn
                disp(errmes)
                disp(error62)
                disp(error6)
            end
            exitflag=-6;
            if trx
                x=x';
            end
            return
        elseif ngc<ZeroGrad
            if dispwarn
                disp(errmes)
                disp(error63)
            end
            exitflag=-6;
            if trx
                x=x';
            end
            return
        end
        g=g+PenCoef*gc; ng=norm(g);
    end
end
grec=g; nng=ng;
% ----}
% INITIAL STEPSIZE
h=h1*sqrt(optim.TolX)*max(abs(x));     % smallest possible stepsize
if abs(optim.minimizer_indicator)~=1
    h=h1*max(abs([optim.minimizer_indicator,h]));     % user-supplied stepsize
else
    h=h1*max(1/log(ng+1.1),abs(h));    % calculated stepsize
end

% RESETTING LOOP ----{
while 1
    kcheck=0;                        % Set checkpoint counter.
    kg=0;                            % stepsizes stored
    kj=0;                            % ravine jump counter
    B=eye(n);                        % re-set transf. matrix to identity
    fst=f; g1=g;  dx=0;
    % ----}

    % MAIN ITERATIONS ----{

    while 1
        k=k+1;kcheck=kcheck+1;
        laststep=dx;

        % ADJUST GAMMA --{
        gamma=1+max([ajb^((ajp-kcheck)*n),2*optim.TolFun]);
        gamma=min([gamma,ajs^max([1,log10(nng+1)])]);
        % --}
        gt=g*B;   w=wdef;
        % JUMPING OVER A RAVINE ----{
        if (gt/norm(gt))*(g1'/norm(g1))<low_bound
            if kj==2
                xx=x;
            end
            if kj==0
                kd=4;
            end
            kj=kj+1;  w=-.9; h=h*2;             % use large coef. of space dilation
            if kj>2*kd
                kd=kd+1;
                warnno=1;
                if any(abs(x-xx)<epsnorm*abs(x)) % flat bottom is detected
                    if dispwarn
                        disp(wrnmes)
                        disp(warn08)
                    end
                end
            end
        else
            kj=0;
        end
        % ----}
        % DILATION ----{
        z=gt-g1;
        nrmz=norm(z);
        if(nrmz>epsnorm*norm(gt))
            z=z/nrmz;
            g1=gt+w*(z*gt')*z;  B=B+w*(B*z')*z;
        else
            z=zeros(1,n);
            nrmz=0;
            g1=gt;
        end
        d1=norm(g1);  g0=(g1/d1)*B';
        % ----}
        % RESETTING ----{
        if kcheck>1
            idx=find(abs(g)>ZeroGrad); numelem=size(idx,2);
            if numelem>0, grbnd=epsnorm*numelem^2;
                if all(abs(g1(idx))<=abs(g(idx))*grbnd) || nrmz==0
                    if dispwarn
                        disp(wrnmes)
                        disp(warn20)
                    end
                    if abs(fst-f)<abs(f)*.01
                        ajp=ajp-10*n;
                    else
                        ajp=ajpp;
                    end
                    h=h1*dx/3;
                    k=k-1;
                    break
                end
            end
        end
        % ----}
        % STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH
        xopt=x;fopt=f;   k1=0;k2=0;ksm=0;kc=0;knan=0;  hp=h;
        if constr, Reset=0; end
        % 1-D SEARCH ----{
        while 1
            x1=x;f1=f;
            if constr
                FP1=FP;
                fp1=fp;
            end
            x=x+hp*g0;
            % FUNCTION VALUE
            if trx
                f=feval(fun,x',varargin{:});
            else
                f=feval(fun,x,varargin{:});
            end
            n_f_evals=n_f_evals+1;
            if h1*f==Inf
                if dispwarn
                    disp(errmes)
                    disp(error5)
                end
                exitflag=-7;
                if trx
                    x=x';
                end
                return
            end
            if constr, fp=f;
                if trx
                    fc=feval(func,x');
                else
                    fc=feval(func,x);
                end
                n_constraint_evals=n_constraint_evals+1;
                if  isnan(fc)
                    if dispwarn
                        disp(errmes)
                        disp(error51)
                        disp(error6)
                    end
                    exitflag=-5;
                    if trx
                        x=x';
                    end
                    return
                elseif abs(fc)==Inf
                    if dispwarn
                        disp(errmes)
                        disp(error52)
                        disp(error6)
                    end
                    exitflag=-5;
                    if trx
                        x=x';
                    end
                    return
                end
                if fc<=cnteps
                    FP=1;
                    fc=0;
                else
                    FP=0;
                    fp_rate=(fp-fp1);
                    if fp_rate<-epsnorm
                        if ~FP1
                            PenCoefNew=-15*fp_rate/norm(x-x1);
                            if PenCoefNew>1.2*PenCoef
                                PenCoef=PenCoefNew; Reset=1; kless=0; f=f+PenCoef*fc; break
                            end
                        end
                    end
                end
                f=f+PenCoef*fc;
            end
            if abs(f)==Inf || isnan(f)
                if dispwarn, disp(wrnmes)
                    if isnan(f)
                        disp(error31)
                    else
                        disp(error32)
                    end
                end
                if ksm || kc>=mxtc
                    exitflag=-3;
                    % don't return with NaN or Inf despite error code
                    x=x1;
                    f=f1;                  
                    if trx
                        x=x';
                    end
                    return
                else
                    k2=k2+1;
                    k1=0;
                    hp=hp/dq;
                    x=x1;
                    f=f1;
                    knan=1;
                    if constr
                        FP=FP1;
                        fp=fp1;
                    end
                end
                % STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
            elseif all(abs(x-x1)<abs(x)*epsnorm)
                stepvanish=stepvanish+1;
                if stepvanish>=5
                    exitflag=-14;
                    if dispwarn
                        disp(termwarn1)
                        disp(endwarn(4,:))
                    end
                    if trx
                        x=x';
                    end
                    return
                else
                    x=x1;
                    f=f1;
                    hp=hp*10;
                    ksm=1;
                    if constr
                        FP=FP1;
                        fp=fp1;
                    end
                end
                % USE SMALLER STEP
            elseif h1*f<h1*gamma^sign(f1)*f1
                if ksm
                    break
                end
                k2=k2+1;k1=0; hp=hp/dq; x=x1;f=f1;
                if constr
                    FP=FP1;
                    fp=fp1;
                end
                if kc>=mxtc, break, end
                % 1-D OPTIMIZER IS LEFT BEHIND
            else
                if h1*f<=h1*f1
                    break
                end
                % USE LARGER STEP
                k1=k1+1;
                if k2>0
                    kc=kc+1;
                end
                k2=0;
                if k1>=20
                    hp=du20*hp;
                elseif k1>=10
                    hp=du10*hp;
                elseif k1>=3
                    hp=du03*hp;
                end
            end
        end
        % ----}  End of 1-D search
        % ADJUST THE TRIAL STEP SIZE ----{
        dx=norm(xopt-x);
        if kg<kstore
            kg=kg+1;
        end
        if kg>=2
            nsteps(2:kg)=nsteps(1:kg-1);
        end
        nsteps(1)=dx/(abs(h)*norm(g0));
        kk=sum(nsteps(1:kg).*[kg:-1:1])/sum(kg:-1:1);
        if     kk>des
            if kg==1
                h=h*(kk-des+1);
            else
                h=h*sqrt(kk-des+1);
            end
        elseif kk<des
            h=h*sqrt(kk/des);
        end

        stepvanish=stepvanish+ksm;
        % ----}
        % COMPUTE THE GRADIENT ----{
        if app
            deltax=sign(g0); idx=find(deltax==0);
            deltax(idx)=ones(size(idx));  deltax=h1*ddx*deltax;
            if constr
                if trx
                    g=apprgrdn(x',fp,fun,deltax',1,varargin{:});
                else
                    g=apprgrdn(x ,fp,fun,deltax,1,varargin{:});
                end
            else
                if trx
                    g=apprgrdn(x',f,fun,deltax',1,varargin{:});
                else
                    g=apprgrdn(x ,f,fun,deltax ,1,varargin{:});
                end
            end
            n_f_evals=n_f_evals+n;
        else
            if trx                
                [~,g]=feval(fun,x',varargin{:});
            else
                [~,g]=feval(fun,x,varargin{:});
            end
            n_grad_evals=n_grad_evals+1;
        end
        if size(g,2)==1
            g=g';
        end
        ng=norm(g);
        if isnan(ng)
            if dispwarn
                disp(errmes)
                disp(error41)
            end
            exitflag=-4;
            if trx
                x=x';
            end
            return
        elseif ng==Inf
            if dispwarn
                disp(errmes)
                disp(error42)
            end
            exitflag=-4;
            if trx
                x=x';
            end
            return
        elseif ng<ZeroGrad
            if dispwarn
                disp(wrnmes)
                disp(warn1)
            end
            ng=ZeroGrad;
        end
        % Constraints:
        if constr
            if ~FP
                if ng<.01*PenCoef
                    kless=kless+1;
                    if kless>=20
                        PenCoef=PenCoef/10;
                        Reset=1;
                        kless=0;
                    end
                else
                    kless=0;
                end
                if appconstr
                    deltax=sign(x); idx=find(deltax==0);
                    deltax(idx)=ones(size(idx));  deltax=ddx*deltax;
                    if trx
                        gc=apprgrdn(x',fc,func,deltax',0);
                    else
                        gc=apprgrdn(x ,fc,func,deltax ,0);
                    end
                    n_constraint_evals=n_constraint_evals+n;
                else
                    if trx
                        gc=feval(gradc,x');
                    else
                        gc=feval(gradc,x );
                    end
                    n_constraint_gradient_evals=n_constraint_gradient_evals+1;
                end
                if size(gc,2)==1
                    gc=gc';
                end
                ngc=norm(gc);
                if isnan(ngc)
                    if dispwarn
                        disp(errmes)
                        disp(error61)
                    end
                    exitflag=-6;
                    if trx
                        x=x';
                    end
                    return
                elseif ngc==Inf
                    if dispwarn
                        disp(errmes)
                        disp(error62)
                    end
                    exitflag=-6;
                    if trx
                        x=x';
                    end
                    return
                elseif ngc<ZeroGrad && ~appconstr
                    if dispwarn
                        disp(errmes)
                        disp(error63)
                    end
                    exitflag=-6;
                    if trx
                        x=x';
                    end
                    return
                end
                g=g+PenCoef*gc; ng=norm(g);
                if Reset
                    if dispwarn
                        disp(wrnmes)
                        disp(warn21)
                    end
                    h=h1*dx/3; k=k-1; nng=ng; break
                end
            end
        end
        if h1*f>h1*frec
            frec=f;
            xrec=x;
            grec=g;
        end
        % ----}
        if ng>ZeroGrad
            if knorms<10
                knorms=knorms+1;
            end
            if knorms>=2
                gnorms(2:knorms)=gnorms(1:knorms-1);
            end
            gnorms(1)=ng;
            nng=(prod(gnorms(1:knorms)))^(1/knorms);
        end

        % DISPLAY THE CURRENT VALUES ----{
        if k==ld
            disp('Iter.# ..... Function ... Step Value ... Gradient Norm ');
            fprintf('%5i   %13.5e   %13.5e     %13.5e\n',k,f,dx,ng);
            ld=k+dispdata;
        end
        %----}
        % CHECK THE STOPPING CRITERIA ----{
        termflag=1;
        if constr
            if ~FP
                termflag=0;
            end
        end
        if kcheck<=5
            termflag=0;
        end
        if knan
            termflag=0;
        end
        if kc>=mxtc
            termflag=0;
        end
        % ARGUMENT
        if termflag
            idx=find(abs(x)>=lowxbound);
            if isempty(idx) || all(abs(xopt(idx)-x(idx))<=optim.TolX*abs(x(idx)))
                termx=termx+1;
                % FUNCTION
                if abs(f-frec)> detfr * abs(f)    && ...
                        abs(f-fopt)<=optim.TolFun*abs(f) && ...
                        krerun<=3                      && ...
                        ~constr
                    if any(abs(xrec(idx)-x(idx))> detxr * abs(x(idx)))
                        if dispwarn
                            disp(wrnmes)
                            disp(warn09)
                        end
                        x=xrec;
                        f=frec;
                        g=grec;
                        ng=norm(g);
                        krerun=krerun+1;
                        h=h1*max([dx,detxr*norm(x)])/krerun;
                        warnno=2;
                        break
                    else
                        h=h*10;
                    end
                elseif  abs(f-frec)> optim.TolFun*abs(f)    && ...
                        norm(x-xrec)<optim.TolX*norm(x) && constr

                elseif  abs(f-fopt)<=optim.TolFun*abs(f)  || ...
                        abs(f)<=lowfbound               || ...
                        (abs(f-fopt)<=optim.TolFun && termx>=limxterm )
                    if stopf
                        if dx<=laststep
                            if warnno==1 && ng<sqrt(optim.TolFun)
                                warnno=0;
                            end
                            if ~app
                                if any(abs(g)<=epsnorm2)
                                    warnno=3;
                                end
                            end
                            if warnno~=0
                                exitflag=-warnno-10;
                                if dispwarn, disp(termwarn1)
                                    disp(endwarn(warnno,:))
                                    if app
                                        disp(appwarn);
                                    end
                                end
                            else
                                exitflag=k;
                                if dispwarn
                                    disp(termwarn0);
                                end
                            end
                            if trx
                                x=x';
                            end
                            return
                        end
                    else
                        stopf=1;
                    end
                elseif dx<1.e-12*max(norm(x),1) && termx>=limxterm
                    exitflag=-14;
                    if dispwarn
                        disp(termwarn1)
                        disp(endwarn(4,:))
                        if app
                            disp(appwarn)
                        end
                    end
                    x=xrec; f=frec;
                    if trx
                        x=x';
                    end
                    return
                else
                    stopf=0;
                end
            end
        end
        % ITERATIONS LIMIT
        if(k==optim.MaxIter)
            exitflag=-9;
            if trx
                x=x';
            end
            if dispwarn
                disp(wrnmes)
                disp(warn4)
            end
            return
        end
        % ----}
        % ZERO GRADIENT ----{
        if constr
            if ng<=ZeroGrad
                if dispwarn
                    disp(termwarn1)
                    disp(warn1)
                end
                exitflag=-8;
                if trx
                    x=x';
                end
                return
            end
        else
            if ng<=ZeroGrad
                nzero=nzero+1;
                if dispwarn
                    disp(wrnmes)
                    disp(warn1)
                end
                if nzero>=3
                    exitflag=-8;
                    if trx
                        x=x';
                    end
                    return
                end
                g0=-h*g0/2;
                for i=1:10
                    x=x+g0;
                    if trx
                        f=feval(fun,x',varargin{:});
                    else
                        f=feval(fun,x,varargin{:});
                    end
                    n_f_evals=n_f_evals+1;
                    if abs(f)==Inf
                        if dispwarn
                            disp(errmes)
                            disp(error32)
                        end
                        exitflag=-3;
                        if trx
                            x=x';
                        end
                        return
                    elseif isnan(f)
                        if dispwarn
                            disp(errmes)
                            disp(error31)
                        end
                        exitflag=-3;
                        if trx
                            x=x';
                        end
                        return
                    end
                    if app
                        deltax=sign(g0);
                        idx=find(deltax==0);
                        deltax(idx)=ones(size(idx));
                        deltax=h1*ddx*deltax;
                        if trx
                            g=apprgrdn(x',f,fun,deltax',1,varargin{:});
                        else
                            g=apprgrdn(x,f,fun,deltax,1,varargin{:});
                        end
                        n_f_evals=n_f_evals+n;
                    else
                        if trx
                            [~,g]=feval(fun,x',varargin{:});
                        else
                            [~,g]=feval(fun,x,varargin{:});
                        end
                        n_grad_evals=n_grad_evals+1;
                    end
                    if size(g,2)==1
                        g=g';
                    end
                    ng=norm(g);
                    if ng==Inf
                        if dispwarn
                            disp(errmes)
                            disp(error42)
                        end
                        exitflag=-4;
                        if trx
                            x=x';
                        end
                        return
                    elseif isnan(ng)
                        if dispwarn
                            disp(errmes)
                            disp(error41)
                        end
                        exitflag=-4;
                        if trx
                            x=x';
                        end
                        return
                    end
                    if ng>ZeroGrad
                        break
                    end
                end
                if ng<=ZeroGrad
                    if dispwarn
                        disp(termwarn1)
                        disp(warn1)
                    end
                    exitflag=-8;
                    if trx
                        x=x';
                    end
                    return
                end
                h=h1*dx;
                break
            end
        end
        % ----}
        % FUNCTION IS FLAT AT THE POINT ----{
        if ~constr && abs(f-fopt)<abs(fopt)*optim.TolFun && kcheck>5 && ng<1
            idx=find(abs(g)<=epsnorm2);
            ni=size(idx,2);
            if ni>=1 && ni<=n/2 && kflat<=3
                kflat=kflat+1;
                if dispwarn
                    disp(wrnmes)
                    disp(warn31)
                end
                warnno=1;
                x1=x; fm=f;
                for j=idx
                    y=x(j); f2=fm;
                    if y==0
                        x1(j)=1;
                    elseif abs(y)<1
                        x1(j)=sign(y);
                    else
                        x1(j)=y;
                    end
                    for i=1:20
                        x1(j)=x1(j)/1.15;
                        if trx
                            f1=feval(fun,x1',varargin{:});
                        else
                            f1=feval(fun,x1,varargin{:});
                        end
                        n_f_evals=n_f_evals+1;
                        if abs(f1)~=Inf && ~isnan(f1)
                            if h1*f1>h1*fm
                                y=x1(j);
                                fm=f1;
                            elseif h1*f2>h1*f1
                                break
                            elseif f2==f1
                                x1(j)=x1(j)/1.5;
                            end
                            f2=f1;
                        end
                    end
                    x1(j)=y;
                end
                if h1*fm>h1*f
                    if app
                        deltax=h1*ddx*ones(size(deltax));
                        if trx
                            gt=apprgrdn(x1',fm,fun,deltax',1,varargin{:});
                        else
                            gt=apprgrdn(x1 ,fm,fun,deltax ,1,varargin{:});
                        end
                        n_f_evals=n_f_evals+n;
                    else
                        if trx
                            [~,gt]=feval(fun,x1',varargin{:});
                        else
                            [~,gt]=feval(fun,x1,varargin{:});
                        end
                        n_grad_evals=n_grad_evals+1;
                    end
                    if size(gt,2)==1
                        gt=gt';
                    end
                    ngt=norm(gt);
                    if ~isnan(ngt) && ngt>epsnorm2
                        if dispwarn
                            disp(warn32)
                        end
                        optim.TolFun=optim.TolFun/5;
                        x=x1;
                        g=gt;
                        ng=ngt;
                        f=fm;
                        h=h1*dx/3;
                        break
                    end
                end
            end
        end
        % ----}
    end % iterations
end % restart
    % end of the function
    %
end