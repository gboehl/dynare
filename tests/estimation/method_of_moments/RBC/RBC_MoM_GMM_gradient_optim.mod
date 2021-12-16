% Test whether gradient-based optimizers are able to use analytical
% Jacobian of moments in GMM estimation
%
% Copyright (C) 2021 Dynare Team
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
  @#endif

  estimated_params_init(use_calibration);
  end;
  
  @#for optimizer in [1, 3, 4, 101, 13]
    @#if estimParams == 2 && optimizer == 13
        %skip due to buggy behavior in Octave
        if ~isoctave
    @#endif
  
  method_of_moments(
          mom_method = GMM         % method of moments method; possible values: GMM|SMM
        , datafile   = 'RBC_Andreasen_Data_2.mat' % name of filename with data
        , order = 2                 % order of Taylor approximation in perturbation
        , weighting_matrix = ['OPTIMAL']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
        , nodisplay
        , nograph
        , mode_compute = @{optimizer}        % specifies the optimizer for minimization of moments distance
        %, additional_optimizer_steps = [1 3 13]
%        , optim = ('DerivativeCheck', 'on','FiniteDifferenceType','central'
%                    ,'TolX', 1e-6
%                    ,'MaxIter', 3000
%                    ,'MaxFunEvals', 1D6
%                    ,'UseParallel' , 1
%                    ,'Jacobian' , 'on'                   
%                    ,'GradObj','on'
%                   )    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
        , silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between        
        , analytic_jacobian
  );
  
    @#if estimParams == 2 && optimizer == 13
        %skip due to buggy behavior in Octave
        end
    @#endif
  @#endfor  
  
@#endfor
