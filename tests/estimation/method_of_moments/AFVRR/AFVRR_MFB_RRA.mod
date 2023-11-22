% DSGE model based on replication files of
% Andreasen, Fernandez-Villaverde, Rubio-Ramirez (2018), The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications, Review of Economic Studies, 85, p. 1-49
% Adapted for Dynare by Willi Mutschler (@wmutschl, willi@mutschler.eu), Jan 2021
% =========================================================================
% Copyright © 2021-2023 Dynare Team
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

% This is the model with feedback and calibrated RRA
% Original code RunGMM_Feedback_estim_RRA_5.m by Martin M. Andreasen, Jan 2016

@#include "AFVRR_common.inc"

%--------------------------------------------------------------------------
% Parameter calibration taken from RunGMM_Feedback_estim_RRA_5.m
%--------------------------------------------------------------------------
% fixed parameters
INHABIT   = 1;
PHI1      = 4;
PHI4      = 1;
KAPAone   = 0;
DELTA     = 0.025;
THETA     = 0.36;
ETA       = 6;
CHI       = 0;
BETTAxhr  = 0;
BETTAxhr40= 0;
RHOD      = 0;
GAMA      = 0.9999;
CONSxhr20 = 0;
RRA       = 5;

% estimated parameters
BETTA     = 0.996850651147000;
B         = 0.684201133923000;
H         = 0.338754441432000;
PHI2      = 0.738293581320000;
KAPAtwo   = 11.664785970704999;
ALFA      = 0.831836572237000;
RHOR      = 0.772754520116000;
BETTAPAI  = 3.020381242896000;
BETTAY    = 0.288367683973000;
MYYPS     = 1.000911709188000;
MYZ       = 1.005433723022000;
RHOA      = 0.749465413198000;
RHOG      = 0.847225569814000;
PAI       = 1.010428794858000;
CONSxhr40 = 0.992863217133000;
GoY       = 0.207099399789000;
STDA      = 0.015621059978000;
STDG      = 0.047539390956000;
STDD      = 0.008623441943000;

% endogenous parameters set via steady state, no need to initialize
%PHIzero  = ;
%AA       = ;
%PHI3     = ;
%negVf    = ;

model_diagnostics;
% Model diagnostics show that some parameters are endogenously determined
% via the steady state, so we run steady to calibrate all parameters
steady;
model_diagnostics;
% Now all parameters are determined

resid;
check;

%--------------------------------------------------------------------------
% Shock distribution
%--------------------------------------------------------------------------
shocks;
var eps_a = STDA^2;
var eps_d = STDD^2;
var eps_g = STDG^2;
end;

%--------------------------------------------------------------------------
% Estimated Params block - these parameters will be estimated, we
% initialize at calibrated values
%--------------------------------------------------------------------------
estimated_params;
BETTA;
B;
H;
PHI2;
KAPAtwo;
ALFA;
RHOR;
BETTAPAI;
BETTAY;
MYYPS;
MYZ;
RHOA;
RHOG;
PAI;
CONSxhr40;
GoY;
stderr eps_a;
stderr eps_g;
stderr eps_d;
end;

estimated_params_init(use_calibration);
end;

%--------------------------------------------------------------------------
% Compare whether toolbox yields equivalent moments at second order
%--------------------------------------------------------------------------
% Note that we compare results for orderApp=1|2 and not for orderApp=3, because
% there is a small error in the replication files of the original article in the 
% computation of the covariance matrix of the extended innovations vector.
% The authors have been contacted, fixed it, and report that the results 
% change only slightly at orderApp=3 to what they report in the paper. At
% orderApp=2 all is correct and so the following part tests whether we get 
% the same model moments at the calibrated parameters (we do not optimize).
% We compare it to the replication file RunGMM_Feedback_estim_RRA.m with the
% following settings: orderApp=1|2, seOn=1, q_lag=10, weighting=1+1;
% scaled=0; optimizer=0; estimator=1; momentSet=2;
%
% Output of the replication files for orderApp=1
AndreasenEtAl.Q1 = 60275.3715;
AndreasenEtAl.moments1 =[ % note that we reshuffeled to be compatible with our matched moments block
    {[ 1]}    {'Ex'  }    {'Gr_C   '}    {'     '  }    {'0.024388'   }    {'0.023726'   }
    {[ 2]}    {'Ex'  }    {'Gr_I   '}    {'     '  }    {'0.031046'   }    {'0.027372'   }
    {[ 3]}    {'Ex'  }    {'Infl  ' }    {'     '  }    {'0.03757'    }    {'0.041499'   }
    {[ 4]}    {'Ex'  }    {'r1    ' }    {'     '  }    {'0.056048'   }    {'0.077843'   }
    {[ 5]}    {'Ex'  }    {'r40   ' }    {'     '  }    {'0.069929'   }    {'0.077843'   }
    {[ 6]}    {'Ex'  }    {'xhr40  '}    {'     '  }    {'0.017237'   }    {'0'          }
    {[ 7]}    {'Ex'  }    {'GoY    '}    {'     '  }    {'-1.5745'    }    {'-1.5746'    }
    {[ 8]}    {'Ex'  }    {'hours  '}    {'     '  }    {'-0.043353'  }    {'-0.043299'  }
    {[ 9]}    {'Exx' }    {'Gr_C   '}    {'Gr_C   '}    {'0.0013159'  }    {'0.0012763'  }
    {[10]}    {'Exx' }    {'Gr_C   '}    {'Gr_I   '}    {'0.0021789'  }    {'0.0017759'  }
    {[11]}    {'Exx' }    {'Gr_C   '}    {'Infl  ' }    {'0.00067495' }    {'0.00077354' }
    {[12]}    {'Exx' }    {'Gr_C   '}    {'r1    ' }    {'0.0011655'  }    {'0.0016538'  }
    {[13]}    {'Exx' }    {'Gr_C   '}    {'r40   ' }    {'0.0015906'  }    {'0.0017949'  }
    {[14]}    {'Exx' }    {'Gr_C   '}    {'xhr40  '}    {'0.0020911'  }    {'0.0017847'  }
    {[15]}    {'Exx' }    {'Gr_I   '}    {'Gr_I   '}    {'0.0089104'  }    {'0.0053424'  }
    {[16]}    {'Exx' }    {'Gr_I   '}    {'Infl  ' }    {'0.00063139' }    {'0.00064897' }
    {[17]}    {'Exx' }    {'Gr_I   '}    {'r1    ' }    {'0.0011031'  }    {'0.0019533'  }
    {[18]}    {'Exx' }    {'Gr_I   '}    {'r40   ' }    {'0.0018445'  }    {'0.0020602'  }
    {[19]}    {'Exx' }    {'Gr_I   '}    {'xhr40  '}    {'0.00095556' }    {'0.0064856'  }
    {[20]}    {'Exx' }    {'Infl  ' }    {'Infl  ' }    {'0.0020268'  }    {'0.0020922'  }
    {[21]}    {'Exx' }    {'Infl  ' }    {'r1    ' }    {'0.0025263'  }    {'0.0036375'  }
    {[22]}    {'Exx' }    {'Infl  ' }    {'r40   ' }    {'0.0029126'  }    {'0.0034139'  }
    {[23]}    {'Exx' }    {'Infl  ' }    {'xhr40  '}    {'-0.00077101'}    {'-0.0011665' }
    {[24]}    {'Exx' }    {'r1    ' }    {'r1    ' }    {'0.0038708'  }    {'0.0066074'  }
    {[25]}    {'Exx' }    {'r1    ' }    {'r40   ' }    {'0.0044773'  }    {'0.0062959'  }
    {[26]}    {'Exx' }    {'r1    ' }    {'xhr40  '}    {'-0.00048202'}    {'-0.00075499'}
    {[27]}    {'Exx' }    {'r40   ' }    {'r40   ' }    {'0.0054664'  }    {'0.0061801'  }
    {[28]}    {'Exx' }    {'r40   ' }    {'xhr40  '}    {'0.00053864' }    {'-0.00030456'}
    {[29]}    {'Exx' }    {'xhr40  '}    {'xhr40  '}    {'0.053097'   }    {'0.012048'   }
    {[30]}    {'Exx' }    {'GoY    '}    {'GoY    '}    {'2.4863'     }    {'2.4872'     }
    {[31]}    {'Exx' }    {'hours  '}    {'hours  '}    {'0.0018799'  }    {'0.0018759'  }
    {[32]}    {'Exx1'}    {'Gr_C   '}    {'Gr_C   '}    {'0.00077917' }    {'0.00080528' }
    {[33]}    {'Exx1'}    {'Gr_I   '}    {'Gr_I   '}    {'0.0050104'  }    {'0.0017036'  }
    {[34]}    {'Exx1'}    {'Infl  ' }    {'Infl  ' }    {'0.0019503'  }    {'0.0020185'  }
    {[35]}    {'Exx1'}    {'r1    ' }    {'r1    ' }    {'0.0038509'  }    {'0.0065788'  }
    {[36]}    {'Exx1'}    {'r40   ' }    {'r40   ' }    {'0.0054699'  }    {'0.0061762'  }
    {[37]}    {'Exx1'}    {'xhr40  '}    {'xhr40  '}    {'-0.00098295'}    {'-4.5519e-15'}
    {[38]}    {'Exx1'}    {'GoY    '}    {'GoY    '}    {'2.4868'     }    {'2.4863'     }
    {[39]}    {'Exx1'}    {'hours  '}    {'hours  '}    {'0.0018799'  }    {'0.0018755'  }
];

% Output of the replication files for orderApp=2
AndreasenEtAl.Q2 = 140.8954;
AndreasenEtAl.moments2 =[ % note that we reshuffeled to be compatible with our matched moments block
    {[ 1]}    {'Ex'  }    {'Gr_C   '}    {'     '  }    {'0.024388'   }    {'0.023726'   }
    {[ 2]}    {'Ex'  }    {'Gr_I   '}    {'     '  }    {'0.031046'   }    {'0.027372'   }
    {[ 3]}    {'Ex'  }    {'Infl  ' }    {'     '  }    {'0.03757'    }    {'0.034618'   }
    {[ 4]}    {'Ex'  }    {'r1    ' }    {'     '  }    {'0.056048'   }    {'0.056437'   }
    {[ 5]}    {'Ex'  }    {'r40   ' }    {'     '  }    {'0.069929'   }    {'0.07051'    }
    {[ 6]}    {'Ex'  }    {'xhr40  '}    {'     '  }    {'0.017237'   }    {'0.014242'   }
    {[ 7]}    {'Ex'  }    {'GoY    '}    {'     '  }    {'-1.5745'    }    {'-1.574'     }
    {[ 8]}    {'Ex'  }    {'hours  '}    {'     '  }    {'-0.043353'  }    {'-0.043351'  }
    {[ 9]}    {'Exx' }    {'Gr_C   '}    {'Gr_C   '}    {'0.0013159'  }    {'0.0012917'  }
    {[10]}    {'Exx' }    {'Gr_C   '}    {'Gr_I   '}    {'0.0021789'  }    {'0.0017862'  }
    {[11]}    {'Exx' }    {'Gr_C   '}    {'Infl  ' }    {'0.00067495' }    {'0.00061078' }
    {[12]}    {'Exx' }    {'Gr_C   '}    {'r1    ' }    {'0.0011655'  }    {'0.0011494'  }
    {[13]}    {'Exx' }    {'Gr_C   '}    {'r40   ' }    {'0.0015906'  }    {'0.0016149'  }
    {[14]}    {'Exx' }    {'Gr_C   '}    {'xhr40  '}    {'0.0020911'  }    {'0.002203'   }
    {[15]}    {'Exx' }    {'Gr_I   '}    {'Gr_I   '}    {'0.0089104'  }    {'0.0054317'  }
    {[16]}    {'Exx' }    {'Gr_I   '}    {'Infl  ' }    {'0.00063139' }    {'0.00045278' }
    {[17]}    {'Exx' }    {'Gr_I   '}    {'r1    ' }    {'0.0011031'  }    {'0.0013672'  }
    {[18]}    {'Exx' }    {'Gr_I   '}    {'r40   ' }    {'0.0018445'  }    {'0.0018557'  }
    {[19]}    {'Exx' }    {'Gr_I   '}    {'xhr40  '}    {'0.00095556' }    {'0.0067742'  }
    {[20]}    {'Exx' }    {'Infl  ' }    {'Infl  ' }    {'0.0020268'  }    {'0.0016583'  }
    {[21]}    {'Exx' }    {'Infl  ' }    {'r1    ' }    {'0.0025263'  }    {'0.0024521'  }
    {[22]}    {'Exx' }    {'Infl  ' }    {'r40   ' }    {'0.0029126'  }    {'0.002705'   }
    {[23]}    {'Exx' }    {'Infl  ' }    {'xhr40  '}    {'-0.00077101'}    {'-0.00065007'}
    {[24]}    {'Exx' }    {'r1    ' }    {'r1    ' }    {'0.0038708'  }    {'0.0038274'  }
    {[25]}    {'Exx' }    {'r1    ' }    {'r40   ' }    {'0.0044773'  }    {'0.004297'   }
    {[26]}    {'Exx' }    {'r1    ' }    {'xhr40  '}    {'-0.00048202'}    {'6.3243e-05' }
    {[27]}    {'Exx' }    {'r40   ' }    {'r40   ' }    {'0.0054664'  }    {'0.0051686'  }
    {[28]}    {'Exx' }    {'r40   ' }    {'xhr40  '}    {'0.00053864' }    {'0.00066645' }
    {[29]}    {'Exx' }    {'xhr40  '}    {'xhr40  '}    {'0.053097'   }    {'0.013543'   }
    {[30]}    {'Exx' }    {'GoY    '}    {'GoY    '}    {'2.4863'     }    {'2.4858'     }
    {[31]}    {'Exx' }    {'hours  '}    {'hours  '}    {'0.0018799'  }    {'0.0018804'  }
    {[32]}    {'Exx1'}    {'Gr_C   '}    {'Gr_C   '}    {'0.00077917' }    {'0.00081772' }
    {[33]}    {'Exx1'}    {'Gr_I   '}    {'Gr_I   '}    {'0.0050104'  }    {'0.0017106'  }
    {[34]}    {'Exx1'}    {'Infl  ' }    {'Infl  ' }    {'0.0019503'  }    {'0.0015835'  }
    {[35]}    {'Exx1'}    {'r1    ' }    {'r1    ' }    {'0.0038509'  }    {'0.0037985'  }
    {[36]}    {'Exx1'}    {'r40   ' }    {'r40   ' }    {'0.0054699'  }    {'0.0051642'  }
    {[37]}    {'Exx1'}    {'xhr40  '}    {'xhr40  '}    {'-0.00098295'}    {'0.00020285' }
    {[38]}    {'Exx1'}    {'GoY    '}    {'GoY    '}    {'2.4868'     }    {'2.4848'     }
    {[39]}    {'Exx1'}    {'hours  '}    {'hours  '}    {'0.0018799'  }    {'0.0018799'  }
];

@#for orderApp in 1:2

method_of_moments(
          mom_method = GMM                    % method of moments method; possible values: GMM|SMM
        , datafile   = 'AFVRR_data.mat'       % name of filename with data
        , bartlett_kernel_lag = 10          % bandwith in optimal weighting matrix
        , order = @{orderApp}                 % order of Taylor approximation in perturbation
        , pruning                             % use pruned state space system at higher-order
        % , verbose                           % display and store intermediate estimation results
        , weighting_matrix = ['DIAGONAL']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename        
        % , TeX                               % print TeX tables and graphics
        % Optimization options that can be set by the user in the mod file, otherwise default values are provided        
        %, huge_number=1D10                   % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
        , mode_compute = 0                    % specifies the optimizer for minimization of moments distance, note that by default there is a new optimizer        
        , optim = ('TolFun', 1e-6
                   ,'TolX', 1e-6
                   ,'MaxIter', 3000
                   ,'MaxFunEvals', 1D6
                   ,'UseParallel' , 1
                   %,'Jacobian' , 'on'
                  )    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
        %, silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
        %, analytic_standard_errors
        , se_tolx=1e-10
);

% Check results

fprintf('****************************************************************\n')
fprintf('Compare Results for perturbation order @{orderApp}\n')
fprintf('****************************************************************\n')
dev_Q            = AndreasenEtAl.Q@{orderApp} - oo_.mom.Q;
dev_datamoments  = str2double(AndreasenEtAl.moments@{orderApp}(:,5)) - oo_.mom.data_moments;
dev_modelmoments = str2double(AndreasenEtAl.moments@{orderApp}(:,6)) - oo_.mom.model_moments;

% There is no table command in Octave
if ~isoctave
table([AndreasenEtAl.Q@{orderApp} ; str2double(AndreasenEtAl.moments@{orderApp}(:,5)) ; str2double(AndreasenEtAl.moments@{orderApp}(:,6))],...
      [oo_.mom.Q                  ; oo_.mom.data_moments                              ; oo_.mom.model_moments                            ],...
      [dev_Q                      ; dev_datamoments                                   ; dev_modelmoments                                 ],...
      'VariableNames', {'Andreasen et al', 'Dynare', 'dev'})
end

if norm(dev_modelmoments)> 1e-4
    warning('Something wrong in the computation of moments at order @{orderApp}')
end

@#endfor

%--------------------------------------------------------------------------
% Replicate estimation at orderApp=3
%--------------------------------------------------------------------------
@#ifdef DoEstimation
method_of_moments(
          mom_method = GMM                    % method of moments method; possible values: GMM|SMM
        , datafile   = 'AFVRR_data.mat'       % name of filename with data
        , bartlett_kernel_lag = 10            % bandwith in optimal weighting matrix
        , order = 3                           % order of Taylor approximation in perturbation
        , pruning                             % use pruned state space system at higher-order
        % , verbose                           % display and store intermediate estimation results
        , weighting_matrix = ['DIAGONAL', 'Optimal']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename        
        % , TeX                               % print TeX tables and graphics
        % Optimization options that can be set by the user in the mod file, otherwise default values are provided        
        %, huge_number=1D10                   % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
        , mode_compute = 13                   % specifies the optimizer for minimization of moments distance, note that by default there is a new optimizer
        , additional_optimizer_steps = [13]
        , optim = ('TolFun', 1e-6
                   ,'TolX', 1e-6
                   ,'MaxIter', 3000
                   ,'MaxFunEvals', 1D6
                   ,'UseParallel' , 1
                   %,'Jacobian' , 'on'
                  )    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
        , silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
        %, analytic_standard_errors
        , se_tolx=1e-10
);
@#endif
