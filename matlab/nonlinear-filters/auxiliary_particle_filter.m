function [LIK,lik] = auxiliary_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% [LIK,lik] = auxiliary_particle_filter(ReducedForm,Y,start,ParticleOptions,ThreadsOptions, options_, M_)
% Evaluates the likelihood of a nonlinear model with the auxiliary particle filter
% allowing eventually resampling.
%
% Copyright © 2011-2023 Dynare Team
%
% This file is part of Dynare (particles module).
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare particles module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% Set default
if isempty(start)
    start = 1;
end
% Get perturbation order
order = options_.order;

% Set flag for prunning
pruning = ParticleOptions.pruning;

% Get steady state and mean.
steadystate = ReducedForm.steadystate;
constant = ReducedForm.constant;
state_variables_steady_state = ReducedForm.state_variables_steady_state;

mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
sample_size = size(Y,2);
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);
number_of_particles = ParticleOptions.number_of_particles;

if ReducedForm.use_k_order_solver
    dr = ReducedForm.dr;
    udr = ReducedForm.udr;
else
    % Set local state space model (first order approximation).
    ghx  = ReducedForm.ghx;
    ghu  = ReducedForm.ghu;
    % Set local state space model (second order approximation).
    ghxx = ReducedForm.ghxx;
    ghuu = ReducedForm.ghuu;
    ghxu = ReducedForm.ghxu;
    ghs2 = ReducedForm.ghs2;
    if (order == 3)
        % Set local state space model (third order approximation).
        ghxxx = ReducedForm.ghxxx;
        ghuuu = ReducedForm.ghuuu;
        ghxxu = ReducedForm.ghxxu;
        ghxuu = ReducedForm.ghxuu;
        ghxss = ReducedForm.ghxss;
        ghuss = ReducedForm.ghuss;
    end
end

% Get covariance matrices
Q = ReducedForm.Q;
H = ReducedForm.H;

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
Q_lower_triangular_cholesky = chol(Q)';

% Set seed for randn().
options_=set_dynare_seed_local_options(options_,'default');

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables+log(det(H));
lik  = NaN(sample_size,1);
LIK  = NaN;

% Initialization of the weights across particles.
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean);
%StateVectors = bsxfun(@plus,zeros(state_variance_rank,number_of_particles),StateVectorMean);
if pruning
    if order == 2
        StateVectors_ = StateVectors;
        state_variables_steady_state_ = state_variables_steady_state;
        mf0_ = mf0;
    elseif order == 3
        StateVectors_ = repmat(StateVectors,3,1);
        state_variables_steady_state_ = repmat(state_variables_steady_state,3,1);
        mf0_ = repmat(mf0,1,3); 
        mask2 = number_of_state_variables+1:2*number_of_state_variables;
        mask3 = 2*number_of_state_variables+1:3*number_of_state_variables;
        mf0_(mask2) = mf0_(mask2)+size(ghx,1);
        mf0_(mask3) = mf0_(mask3)+2*size(ghx,1);
    else
        error('Pruning is not available for orders > 3');
    end
end

for t=1:sample_size
    yhat = bsxfun(@minus,StateVectors,state_variables_steady_state);
    if pruning
        yhat_ = bsxfun(@minus,StateVectors_,state_variables_steady_state_);
        if order == 2
            tmp = local_state_space_iteration_2(yhat,zeros(number_of_structural_innovations,number_of_particles),ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,steadystate,ThreadsOptions.local_state_space_iteration_2);
        elseif order == 3
            tmp = local_state_space_iteration_3(yhat_, zeros(number_of_structural_innovations,number_of_particles), ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, pruning);
        else
            error('Pruning is not available for orders > 3');
        end
    else
        if ReducedForm.use_k_order_solver
            tmp = local_state_space_iteration_k(yhat, zeros(number_of_structural_innovations,number_of_particles), dr, M_, options_, udr);
        else
            if order == 2
                tmp = local_state_space_iteration_2(yhat,zeros(number_of_structural_innovations,number_of_particles),ghx,ghu,constant,ghxx,ghuu,ghxu,ThreadsOptions.local_state_space_iteration_2);
            elseif order == 3
                tmp = local_state_space_iteration_3(yhat, zeros(number_of_structural_innovations,number_of_particles), ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, pruning);
            else
                error('Order > 3: use_k_order_solver should be set to true');
            end
        end
    end
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    z = sum(PredictionError.*(H\PredictionError),1) ;
    tau_tilde = weights.*(tpdf(z,3*ones(size(z)))+1e-99) ;
    tau_tilde = tau_tilde/sum(tau_tilde) ;
    indx = resample(0,tau_tilde',ParticleOptions);
    if pruning
        yhat_ = yhat_(:,indx) ;
    end
    yhat = yhat(:,indx) ;
    weights_stage_1 = weights(indx)./tau_tilde(indx) ;
    epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    if pruning
        if order == 2
            [tmp, tmp_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,steadystate,ThreadsOptions.local_state_space_iteration_2);
        elseif order == 3
            [tmp, tmp_] = local_state_space_iteration_3(yhat_, epsilon, ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, pruning);
        else
            error('Pruning is not available for orders > 3');
        end
        StateVectors_ = tmp_(mf0_,:);
    else
        if ReducedForm.use_k_order_solver
            tmp = local_state_space_iteration_k(yhat, epsilon, dr, M_, options_, udr);
        else
            if order == 2
                tmp = local_state_space_iteration_2(yhat, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, ThreadsOptions.local_state_space_iteration_2);
            elseif order == 3
                tmp = local_state_space_iteration_3(yhat, epsilon, ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, pruning);
            else
                error('Order > 3: use_k_order_solver should be set to true');
            end
        end
    end
    StateVectors = tmp(mf0,:);
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    weights_stage_2 = weights_stage_1.*(exp(-.5*(const_lik+sum(PredictionError.*(H\PredictionError),1))) + 1e-99) ;
    lik(t) = log(mean(weights_stage_2)) ;
    weights = weights_stage_2/sum(weights_stage_2);
    if (ParticleOptions.resampling.status.generic && neff(weights)<ParticleOptions.resampling.threshold*sample_size) || ParticleOptions.resampling.status.systematic
        if pruning
            temp = resample([StateVectors' StateVectors_'],weights',ParticleOptions);
            StateVectors = temp(:,1:number_of_state_variables)';
            StateVectors_ = temp(:,number_of_state_variables+1:end)';
        else
            StateVectors = resample(StateVectors',weights',ParticleOptions)';
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    end
end

LIK = -sum(lik(start:end));
