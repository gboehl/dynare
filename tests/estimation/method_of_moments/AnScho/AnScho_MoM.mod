% DSGE model used in replication files of
% An, Sungbae and Schorfheide, Frank, (2007), Bayesian Analysis of DSGE Models, Econometric Reviews, 26, issue 2-4, p. 113-172.
% Adapted by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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

% Define testscenario
@#define orderApp = 2
@#define estimParams = 1

% Note that we set the numerical optimization tolerance levels very large to speed up the testsuite
@#define optimizer = 13

var c p R g y z INFL INT YGR;
varexo e_r e_g e_z;
parameters tau nu kap cyst psi1 psi2 rhor rhog rhoz rrst pist gamst;

varobs INT YGR INFL;

tau   = 2;
nu    = 0.1;
kap   = 0.33;
cyst  = 0.85;
psi1  = 1.5;
psi2  = 0.125;
rhor  = 0.75;
rhog  = 0.95;
rhoz  = 0.9;
rrst  = 1;
pist  = 3.2;
gamst = 0.55;

model;
#pist2 = exp(pist/400);
#rrst2 = exp(rrst/400);
#bet   = 1/rrst2;
#phi   = tau*(1-nu)/nu/kap/pist2^2;
#gst   = 1/cyst;
#cst   = (1-nu)^(1/tau);
#yst   = cst*gst;
#dy    = y-y(-1);
1 = exp(-tau*c(+1)+tau*c+R-z(+1)-p(+1));
(1-nu)/nu/phi/(pist2^2)*(exp(tau*c)-1) = (exp(p)-1)*((1-1/2/nu)*exp(p)+1/2/nu) - bet*(exp(p(+1))-1)*exp(-tau*c(+1)+tau*c+y(+1)-y+p(+1));
exp(c-y) = exp(-g) - phi*pist2^2*gst/2*(exp(p)-1)^2;
R = rhor*R(-1) + (1-rhor)*psi1*p + (1-rhor)*psi2*(dy+z) + e_r/100;
g = rhog*g(-1) + e_g/100;
z = rhoz*z(-1) + e_z/100;
YGR = gamst+100*(dy+z);
INFL = pist+400*p;
INT = pist+rrst+4*gamst+400*R;
end;

steady_state_model;
  z = 0; p = 0; g = 0; r = 0; c = 0; y = 0;
  YGR = gamst; INFL = pist; INT = pist + rrst + 4*gamst;
end;

shocks;
  var e_r = 0.20^2;
  var e_g = 0.80^2;
  var e_z = 0.45^2;
  corr e_r,e_g = 0.2;
end;

@#if estimParams == 0
% Define only initial values without bounds
estimated_params;
  %tau,             1.50;
  %kap,             0.15;
  psi1,            1.20;
  psi2,            0.50;
  rhor,            0.50;
  %rhog,            0.50;
  %rhoz,            0.50;
  %rrst,            1.20;
  %pist,            3.00;
  gamst,           0.75;
  stderr e_r,      0.30;
  stderr e_g,      0.30;
  stderr e_z,      0.30;
  corr e_r,e_g,    0.10;
end;
@#endif

@#if estimParams == 1
% Define initial values and bounds
estimated_params;
  %tau,             1.50,         1e-5,        10;
  %kap,             0.15,         1e-5,        10;
  psi1,            1.20,         1e-5,        10;
  psi2,            0.50,         1e-5,        10;
  rhor,            0.50,         1e-5,        0.99999;
  %rhog,            0.50,         1e-5,        0.99999;
  %rhoz,            0.50,         1e-5,        0.99999;
  %rrst,            1.20,         1e-5,        10;
  %pist,            3.00,         1e-5,        20;
  gamst,           0.75,         -5,          5;
  stderr e_r,      0.30,         1e-8,        5;
  stderr e_g,      0.30,         1e-8,        5;
  stderr e_z,      0.30,         1e-8,        5;
  corr e_r,e_g,    0.10,         -1,          1;
end;
@#endif

@#if estimParams == 2
% Define prior distribution
estimated_params;
  %tau,             1.50,          1e-5,        10,          gamma_pdf,     2.00,       0.50;
  %kap,             0.15,          1e-5,        10,          gamma_pdf,     0.33,       0.10;
  psi1,            1.20,          1e-5,        10,          gamma_pdf,     1.50,       0.25;
  psi2,            0.50,          1e-5,        10,          gamma_pdf,     0.125,      0.25;
  rhor,            0.50,          1e-5,        0.99999,     beta_pdf,      0.50,       0.20;
  %rhog,            0.50,          1e-5,        0.99999,     beta_pdf,      0.80,       0.10;
  %rhoz,            0.50,          1e-5,        0.99999,     beta_pdf,      0.66,       0.15;
  %rrst,            1.20,          1e-5,        10,          gamma_pdf,     0.50,       0.50;
  %pist,            3.00,          1e-5,        20,          gamma_pdf,     7.00,       2.00;
  gamst,           0.75,          -5,          5,           normal_pdf,    0.40,       0.20;
  stderr e_r,      0.30,          1e-8,        5,           inv_gamma_pdf, 0.50,       0.26;
  stderr e_g,      0.30,          1e-8,        5,           inv_gamma_pdf, 1.25,       0.65;
  stderr e_z,      0.30,          1e-8,        5,           inv_gamma_pdf, 0.63,       0.33;
  corr e_r,e_g,    0.10,          -1,          1,           uniform_pdf,       ,           , -1, 1;
end;
@#endif


% Simulate data
stoch_simul(order=@{orderApp},pruning,nodisplay,nomoments,periods=750,drop=500);
save('AnScho_MoM_data_@{orderApp}.mat', options_.varobs{:} );
pause(1);


%--------------------------------------------------------------------------
% Method of Moments Estimation
%--------------------------------------------------------------------------
matched_moments;
YGR;
INFL;
INT;
%second-order contemporenous product moments
YGR*YGR;
YGR*INFL;
YGR*INT;
INFL*INFL;
INFL*INT;
INT*INT;
%second-order temporal product moments
YGR*YGR(-1);
INT*INT(-1);
INFL*INFL(-1);
end;

% get indices in declaration order
iYGR  = strmatch('YGR',  M_.endo_names,'exact');
iINFL = strmatch('INFL', M_.endo_names,'exact');
iINT  = strmatch('INT',  M_.endo_names,'exact');
% first entry: number of variable in declaration order
% second entry: lag
% third entry: power

matched_moments_ = {
    %first-order product moments
    [iYGR       ]  [0   ],  [1  ];
    [iINFL      ]  [0   ],  [1  ];
    [iINT       ]  [0   ],  [1  ];
    %second-order contemporenous product moments
    [iYGR  iYGR ]  [0  0],  [1 1];
    [iYGR  iINFL]  [0  0],  [1 1];
    [iYGR  iINT ]  [0  0],  [1 1];
    [iINFL iINFL]  [0  0],  [1 1];
    [iINFL iINT ]  [0  0],  [1 1];
    [iINT  iINT ]  [0  0],  [1 1];
    %second-order temporal product moments
    [iYGR  iYGR ]  [0 -1],  [1 1];
    %[iINT  iYGR ]  [0 -1],  [1 1];
    %[iINFL iYGR ]  [0 -1],  [1 1];
    %[iYGR  iINT ]  [0 -1],  [1 1];
    [iINT  iINT ]  [0 -1],  [1 1];
    %[iINFL iINT ]  [0 -1],  [1 1];
    %[iYGR  iINFL]  [0 -1],  [1 1];
    %[iINT  iINFL]  [0 -1],  [1 1];
    [iINFL iINFL]  [0 -1],  [1 1];
};

if ~isequal(M_.matched_moments,matched_moments_)
    error('Translation to matched_moments-block failed')
end

@#for mommethod in ["GMM", "SMM"]
    method_of_moments(
    % Necessery options
          mom_method = @{mommethod}           % method of moments method; possible values: GMM|SMM
        , datafile   = 'AnScho_MoM_data_@{orderApp}.mat' % name of filename with data

    % Options for both GMM and SMM
        % , bartlett_kernel_lag = 20          % bandwith in optimal weighting matrix
        , order = @{orderApp}                 % order of Taylor approximation in perturbation
        % , penalized_estimator               % include deviation from prior mean as additional moment restriction and use prior precision as weight
        , pruning                             % use pruned state space system at higher-order
        % , verbose                           % display and store intermediate estimation results
        , weighting_matrix = ['optimal']      % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
        %, weighting_matrix_scaling_factor=1  % scaling of weighting matrix in objective function
        , se_tolx=1e-6                        % step size for numerical computation of standard errors

    % Options for SMM
        % , burnin=500                        % number of periods dropped at beginning of simulation
        % , bounded_shock_support             % trim shocks in simulation to +- 2 stdev
        % , seed = 24051986                   % seed used in simulations
        % , simulation_multiple = 5           % multiple of the data length used for simulation

    % Options for GMM
        @#if mommethod == "GMM"
        , analytic_standard_errors            % compute standard errors using analytical derivatives
        @#endif

    % General options
        % , dirname = 'MM'                    % directory in which to store estimation output
        % , graph_format = EPS                % specify the file format(s) for graphs saved to disk
        % , nodisplay                         % do not display the graphs, but still save them to disk
        % , nograph                           % do not create graphs (which implies that they are not saved to the disk nor displayed)
        % , noprint                           % do not print stuff to console
        % , plot_priors = 1                   % control plotting of priors
        % , prior_trunc = 1e-10               % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
        % , TeX                               % print TeX tables and graphics

    % Data and model options
        % , first_obs = 501                   % number of first observation
        % , logdata                           % if data is already in logs
        , nobs = 250                          % number of observations
        % , prefilter=0                       % demean each data series by its empirical mean and use centered moments

        % , xls_sheet = data                  % name/number of sheet with data in Excel
        % , xls_range = B2:D200               % range of data in Excel sheet

    % Optimization options that can be set by the user in the mod file, otherwise default values are provided
        % , huge_number=1e7                   % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
        , mode_compute = @{optimizer}         % specifies the optimizer for minimization of moments distance
        , additional_optimizer_steps = [1]    % vector of additional mode-finders run after mode_compute
        % optim: a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute, some exemplary common options:
        , optim = ('TolFun'      , 1e-6       % termination tolerance on the function value, a positive scalar
                  ,'TolX'        , 1e-6       % termination tolerance on x, a positive scalar
                  ,'MaxIter'     , 3000       % maximum number of iterations allowed, a positive integer
                  ,'MaxFunEvals' , 1D6        % maximum number of function evaluations allowed, a positive integer
        %           ,'UseParallel' , 1        % when true (and supported by optimizer) solver estimates gradients in parallel (using Matlab/Octave's parallel toolbox)
        %           ,'Jacobian'    , 'off'    % when 'off' gradient-based solvers approximate Jacobian using finite differences; for GMM we can also pass the analytical Jacobian to gradient-based solvers by setting this 'on'
                  )                         
        , silent_optimizer                    % run minimization of moments distance silently without displaying results or saving files in between

    % Numerical algorithms options
        % , aim_solver                             % Use AIM algorithm to compute perturbation approximation
        % , k_order_solver                         % use k_order_solver in higher order perturbation approximations
        % , dr=default                             % method used to compute the decision rule; possible values are DEFAULT, CYCLE_REDUCTION, LOGARITHMIC_REDUCTION
        % , dr_cycle_reduction_tol = 1e-7          % convergence criterion used in the cycle reduction algorithm
        % , dr_logarithmic_reduction_tol = 1e-12   % convergence criterion used in the logarithmic reduction algorithm
        % , dr_logarithmic_reduction_maxiter = 100 % maximum number of iterations used in the logarithmic reduction algorithm
        % , lyapunov = DEFAULT                     % algorithm used to solve lyapunov equations; possible values are DEFAULT, FIXED_POINT, DOUBLING, SQUARE_ROOT_SOLVER
        % , lyapunov_complex_threshold = 1e-15     % complex block threshold for the upper triangular matrix in symmetric Lyapunov equation solver
        % , lyapunov_fixed_point_tol = 1e-10       % convergence criterion used in the fixed point Lyapunov solver
        % , lyapunov_doubling_tol = 1e-16          % convergence criterion used in the doubling algorithm
        % , sylvester = default                    % algorithm to solve Sylvester equation; possible values are DEFAULT, FIXED_POINT
        % , sylvester_fixed_point_tol = 1e-12      % convergence criterion used in the fixed point Sylvester solver
        % , qz_criterium = 0.999999                % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems
        % , qz_zero_threshold = 1e-6               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition
        % , schur_vec_tol=1e-11                    % tolerance level used to find nonstationary variables in Schur decomposition of the transition matrix
        % , mode_check                             % plot the target function for values around the computed minimum for each estimated parameter in turn
        % , mode_check_neighbourhood_size = 5      % width of the window (expressed in percentage deviation) around the computed minimum to be displayed on the diagnostic plots
        % , mode_check_symmetric_plots=1           % ensure that the check plots are symmetric around the minimum
        % , mode_check_number_of_points = 20       % number of points around the minimum where the target function is evaluated (for each parameter)
    );
@#endfor

