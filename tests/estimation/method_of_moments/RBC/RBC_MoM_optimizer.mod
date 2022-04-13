% Test optimizers
%
% Copyright © 2020-2021 Dynare Team
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
% =========================================================================
% TO DO
% [ ] fix optimizers 11 and 12; 
% note that 12 and 102 require GADS_Toolbox which is not available on servers, but need to be tested locally

% Define testscenario
@#define orderApp = 2

% Note that we will set the numerical optimization tolerance levels very large to speed up the testsuite

@#include "RBC_MoM_common.inc"

shocks;
var u_a; stderr 0.0072;        
end;

varobs c iv n;

%--------------------------------------------------------------------------
% Method of Moments Estimation
%--------------------------------------------------------------------------
matched_moments;
c;
n;
iv;
c*c;
c*iv;
iv*n;
iv*iv;
n*c;
n*n;
c*c(-1);
n*n(-1);
iv*iv(-1);
end;

% reduce options to speed up testsuite
options_.newrat.maxiter = 10;
options_.newrat.tolerance.f = 1e-2;
options_.newrat.tolerance.f_analytic = 1e-2;

options_.mh_jscale = 0.6;
options_.gmhmaxlik.iterations=1;
options_.gmhmaxlik.number=2000;
options_.gmhmaxlik.nclimb=2000;
options_.gmhmaxlik.nscale=2000;
options_.gmhmaxlik.target=0.5;

options_.solveopt.MaxIter=300;
options_.solveopt.LBGradientStep=1e-3;
options_.solveopt.TolFun = 1e-3;
options_.solveopt.TolX = 1e-3;
options_.solveopt.TolXConstraint=1e-3;


@#for estimParams in [0, 1, 2]
  @#if estimParams == 0
    estimated_params(overwrite);
        %DELTA,         0.025;
        %BETTA,         0.984;
        %B,             0.5;
        %ETAc,          2;
        ALFA,          0.667;
        RHOA,          0.979;
        stderr u_a,    0.0072;
    end;
  @#define OPTIMIZERS = [1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 101]
  @#endif

  @#if estimParams == 1
    estimated_params(overwrite);
        %DELTA,         ,        0,           1;
        %BETTA,         ,        0,           1;
        %B,             ,        0,           1;
        %ETAc,          ,        0,           10;
        ALFA,          ,        0,           1;
        RHOA,          ,        0,           1;
        stderr u_a,    ,        0,           1;
    end;
  @#define OPTIMIZERS = [1, 2, 3, 4, 7, 8, 9, 10, 13, 101]
  @#endif

  @#if estimParams == 2
    estimated_params(overwrite);
        %DELTA,         0.025,         0,           1,  normal_pdf, 0.02, 0.5;
        %BETTA,         0.98,         0,           1,  beta_pdf, 0.90, 0.25;
        %B,             0.45,         0,           1,  normal_pdf, 0.40, 0.5;
        %ETAl,          1,            0,           10, normal_pdf, 0.25, 0.0.1;
        %ETAc,          1.8,         0,           10, normal_pdf, 1.80, 0.5;
        ALFA,          0.65,         0,           1,  normal_pdf, 0.60, 0.5;
        RHOA,          0.95,         0,           1,  normal_pdf, 0.90, 0.5;
        stderr u_a,    0.01,         0,           1,  normal_pdf, 0.01, 0.5;
        %THETA,         3.48,          0,           10, normal_pdf, 0.25, 0.0.1;
    end;
  @#define OPTIMIZERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 101]
  @#endif

    estimated_params_init(use_calibration);
    end;

    @#for optimizer in OPTIMIZERS
        
        @#if estimParams == 2 && optimizer == 13
            %skip due to buggy behavior in Octave
        if ~isoctave
        @#endif
        
    method_of_moments(
          mom_method = GMM         % method of moments method; possible values: GMM|SMM
        , datafile   = 'RBC_Andreasen_Data_2.mat' % name of filename with data
        , order = @{orderApp}                 % order of Taylor approximation in perturbation
        , weighting_matrix = ['OPTIMAL']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
        , nograph                           % do not create graphs (which implies that they are not saved to the disk nor displayed)
        , mode_compute = @{optimizer}        % specifies the optimizer for minimization of moments distance
        @#if optimizer == 102
        , optim = ('TolFun'      , 1D-3       % termination tolerance on the function value, a positive scalar
                  ,'MaxIter'     , 300        % maximum number of iterations allowed, a positive integer
                  ,'MaxFunEvals' , 1D3        % maximum number of function evaluations allowed, a positive integer
                  )
        @#else
        , optim = ('TolFun'      , 1D-3       % termination tolerance on the function value, a positive scalar
                  ,'TolX'        , 1e-3       % termination tolerance on x, a positive scalar
                  ,'MaxIter'     , 300        % maximum number of iterations allowed, a positive integer
                  ,'MaxFunEvals' , 1D3        % maximum number of function evaluations allowed, a positive integer
                  )
        @#endif
        %, silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
    );
    
        @#if estimParams == 2 && optimizer == 13
            %skip due to buggy behavior in Octave
        end
        @#endif
    @#endfor
@#endfor



