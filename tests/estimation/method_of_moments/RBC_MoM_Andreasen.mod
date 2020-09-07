% Tests SMM and GMM routines
%
% Copyright (C) 2020 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

% Define testscenario
@#define orderApp = 2
@#define estimParams = 1

% Note that we will set the numerical optimization tolerance levels very large to speed up the testsuite
@#define optimizer = 13


@#include "RBC_MoM_common.inc"

shocks;
var u_a; stderr 0.0072;        
end;

varobs c iv n;


@#if estimParams == 0
estimated_params;
    DELTA,         0.025;
    BETTA,         0.984;
    B,             0.5;
    ETAc,          2;
    ALFA,          0.667;
    RHOA,          0.979;
    stderr u_a,    0.0072;
end;
@#endif

@#if estimParams == 1
estimated_params;
    DELTA,         ,        0,           1;
    BETTA,         ,        0,           1;
    B,             ,        0,           1;
    ETAc,          ,        0,           10;
    ALFA,          ,        0,           1;
    RHOA,          ,        0,           1;
    stderr u_a,    ,        0,           1;
end;
@#endif

@#if estimParams == 2
estimated_params;
    DELTA,         0.025,         0,           1,  normal_pdf, 0.02, 0.5;
    BETTA,         0.98,         0,           1,  beta_pdf, 0.90, 0.25;
    B,             0.45,         0,           1,  normal_pdf, 0.40, 0.5;
    %ETAl,          1,            0,           10, normal_pdf, 0.25, 0.0.1;
    ETAc,          1.8,         0,           10, normal_pdf, 1.80, 0.5;
    ALFA,          0.65,         0,           1,  normal_pdf, 0.60, 0.5;
    RHOA,          0.95,         0,           1,  normal_pdf, 0.90, 0.5;
    stderr u_a,    0.01,         0,           1,  normal_pdf, 0.01, 0.5;
    %THETA,         3.48,          0,           10, normal_pdf, 0.25, 0.0.1;
end;
@#endif

% Simulate data
%stoch_simul(order=@{orderApp},pruning,nodisplay,nomoments,periods=500);
%save('RBC_MoM_data_@{orderApp}.mat', options_.varobs{:} );
%pause(1);


estimated_params_init(use_calibration);
end;


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

c*c(-3);
n*n(-3);
iv*iv(-3);

c*c(-5);
n*n(-5);
iv*iv(-5);
end;

% get indices in declaration order
ic  = strmatch('c',  M_.endo_names,'exact');
iiv = strmatch('iv', M_.endo_names,'exact');
in  = strmatch('n',  M_.endo_names,'exact');
% first entry: number of variable in declaration order
% second entry: lag
% third entry: power

matched_moments_ = {
    [ic     ]  [0   ],  [1  ];
    [in     ]  [0   ],  [1  ];    
    [iiv    ]  [0   ],  [1  ];
    
    [ic  ic ]  [0  0],  [1 1];
    [ic  iiv]  [0  0],  [1 1];
    %[ic  in ]  [0  0],  [1 1];
    %[iiv ic ]  [0  0],  [1 1];
    [in iiv]  [0  0],  [1 1];
    [iiv iiv]  [0  0],  [1 1];    
    [ic  in]  [0  0],  [1 1];
    %[in  iiv]  [0  0],  [1 1];
    [in  in ]  [0  0],  [1 1];
    
    [ic  ic ]  [0 -1],  [1 1];
    [in  in ]  [0 -1],  [1 1];
    [iiv iiv]  [0 -1],  [1 1];

    [ic  ic ]  [0 -3],  [1 1];
    [in  in ]  [0 -3],  [1 1];
    [iiv iiv]  [0 -3],  [1 1];

    [ic  ic ]  [0 -5],  [1 1];
    [in  in ]  [0 -5],  [1 1];
    [iiv iiv]  [0 -5],  [1 1];

};

if ~isequal(M_.matched_moments,matched_moments_)
    error('Translation to matched_moments-block failed')
end


    method_of_moments(
        % Necessery options
          mom_method = GMM                  % method of moments method; possible values: GMM|SMM
        , datafile   = 'RBC_Andreasen_Data_2.mat'         % name of filename with data

        % Options for both GMM and SMM
        %, bartlett_kernel_lag = 20          % bandwith in optimal weighting matrix
        , order = @{orderApp}                % order of Taylor approximation in perturbation
        %, penalized_estimator               % use penalized optimization
        %, pruning                             % use pruned state space system at higher-order
        %, verbose                           % display and store intermediate estimation results
        , weighting_matrix = ['DIAGONAL','OPTIMAL']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename
        %, weighting_matrix_scaling_factor=1
        , additional_optimizer_steps = [13]    % vector of additional mode-finders run after mode_compute
        %, prefilter=0                       % demean each data series by its empirical mean and use centered moments
        % 
        % Options for SMM
        %, bounded_shock_support             % trim shocks in simulation to +- 2 stdev
        %, drop = 500                        % number of periods dropped at beginning of simulation
        %, seed = 24051986                   % seed used in simulations
        %, simulation_multiple = 5           % multiple of the data length used for simulation
        %, burnin = 200
        % 
        % General options
        %, dirname = 'MM'                    % directory in which to store estimation output
        %, graph_format = EPS                % specify the file format(s) for graphs saved to disk
        %, nodisplay                         % do not display the graphs, but still save them to disk
        %, nograph                           % do not create graphs (which implies that they are not saved to the disk nor displayed)
        %, noprint                           % do not print stuff to console
        %, plot_priors = 1                   % control plotting of priors
        %, prior_trunc = 1e-10               % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
        , TeX                               % print TeX tables and graphics
        % 
        % Data and model options
        %, first_obs = 501                     % number of first observation
        %, logdata                           % if loglinear is set, this option is necessary if the user provides data already in logs, otherwise the log transformation will be applied twice (this may result in complex data)
        %, loglinear                         % computes a log-linear approximation of the model instead of a linear approximation
        %, nobs = 50                        % number of observations
        % , xls_sheet = willi                 % name of sheet with data in Excel
        % , xls_range = B2:D200               % range of data in Excel sheet
        % 
        % Optimization options that can be set by the user in the mod file, otherwise default values are provided
        %, analytic_derivation               % uses analytic derivatives to compute standard errors for GMM
        %, huge_number=1D10                   % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
        , mode_compute = 13         % specifies the optimizer for minimization of moments distance, note that by default there is a new optimizer
        , optim = ('TolFun', 1D-6
                   ,'TolX', 1D-6
                  )    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
        %, silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
        , se_tolx = 1e-6                       % convergence criterion on funciton input for numerical differentiation
        % 
        % % Numerical algorithms options
        %, aim_solver                             % Use AIM algorithm to compute perturbation approximation
        %, dr=DEFAULT                             % method used to compute the decision rule; possible values are DEFAULT, CYCLE_REDUCTION, LOGARITHMIC_REDUCTION
        %, dr_cycle_reduction_tol = 1e-7          % convergence criterion used in the cycle reduction algorithm
        %, dr_logarithmic_reduction_maxiter = 100 % maximum number of iterations used in the logarithmic reduction algorithm
        %, dr_logarithmic_reduction_tol = 1e-12   % convergence criterion used in the cycle reduction algorithm
        %, k_order_solver                         % use k_order_solver in higher order perturbation approximations
        %, lyapunov = DEFAULT                     % algorithm used to solve lyapunov equations; possible values are DEFAULT, FIXED_POINT, DOUBLING, SQUARE_ROOT_SOLVER
        %, lyapunov_complex_threshold = 1e-15     % complex block threshold for the upper triangular matrix in symmetric Lyapunov equation solver
        %, lyapunov_fixed_point_tol = 1e-10       % convergence criterion used in the fixed point Lyapunov solver
        %, lyapunov_doubling_tol = 1e-16          % convergence criterion used in the doubling algorithm
        %, sylvester = default                    % algorithm to solve Sylvester equation; possible values are DEFAULT, FIXED_POINT
        %, sylvester_fixed_point_tol = 1e-12      % convergence criterion used in the fixed point Sylvester solver
        %, qz_criterium = 0.999999                % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems [IS THIS CORRET @wmutschl]
        %, qz_zero_threshold = 1e-6               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition
        , mode_check
        %, mode_check_neighbourhood_size=0.5
        %, mode_check_symmetric_plots=0
        %, mode_check_number_of_points=25
    );


