function mdd = hssmc(TargetFun, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_)

% Sequential Monte-Carlo sampler, Herbst and Schorfheide (JAE, 2014).
%
% INPUTS
% - TargetFun        [char]     string specifying the name of the objective function (posterior kernel).
% - xparam1          [double]   p×1 vector of parameters to be estimated (initial values).
% - mh_bounds        [double]   p×2 matrix defining lower and upper bounds for the parameters.
% - dataset_         [dseries]  sample
% - dataset_info     [struct]   informations about the dataset
% - options_         [struct]   dynare's options
% - M_               [struct]   model description
% - estim_params_    [struct]   estimated parameters
% - bayestopt_       [struct]   estimated parameters
% - oo_              [struct]   outputs
%
% SPECIAL REQUIREMENTS
% None.

% Copyright © 2022-2023 Dynare Team
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

    smcopt = options_.posterior_sampler_options.current_options;

    % Set location for the simulated particles.
    SimulationFolder = CheckPath('hssmc', M_.dname);

    % Define prior distribution
    Prior = dprior(bayestopt_, options_.prior_trunc);

    % Set function handle for the objective
    eval(sprintf('%s = @(x) %s(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, mh_bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state, []);', 'funobj', func2str(TargetFun)));

    mlogit = @(x) .95 + .1/(1 + exp(-16*x)); % Update of the scale parameter

    % Create the tempering schedule
    phi = ((0:smcopt.steps-1)/(smcopt.steps-1)).^smcopt.lambda;

    % Initialise the estimate of the marginal density of the data
    mdd = .0;

    % tuning for MH algorithms matrices
    scl = zeros(smcopt.steps, 1);        % scale parameter
    ESS = zeros(smcopt.steps, 1);        % ESS
    acpt = zeros(smcopt.steps, 1);       % average acceptance rate

    % Initialization of the sampler (draws from the prior distribution with finite logged likelihood)
    t0 = tic;
    [particles, tlogpostkernel, loglikelihood] = ...
        smc_samplers_initialization(funobj, 'hssmc', smcopt.particles, Prior, SimulationFolder, smcopt.steps);
    tt = toc(t0);

    dprintf('#Iter.       lambda        ESS          Acceptance rate   scale   resample   seconds')
    dprintf('%3u          %5.4f     %9.5E         %5.4f        %4.3f     %3s       %5.2f', 1, 0, 0, 0, 0, 'no', tt)

    weights = ones(smcopt.particles, 1)/smcopt.particles;

    resampled_particle_swarm = false;

    for i=2:smcopt.steps % Loop over the weight on the liklihood (phi)
        weights = weights.*exp((phi(i)-phi(i-1))*loglikelihood);
        sweight = sum(weights);
        weights = weights/sweight;
        mdd = mdd + log(sweight);
        ESS(i) = 1.0/sum(weights.^2);
        if (2*ESS(i) < smcopt.particles) % Resampling
            resampled_particle_swarm = true;
            iresample = kitagawa(weights);
            particles = particles(:,iresample);
            loglikelihood = loglikelihood(iresample);
            tlogpostkernel = tlogpostkernel(iresample);
            weights = ones(smcopt.particles, 1)/smcopt.particles;
        end
        smcopt.scale = smcopt.scale*mlogit(smcopt.acpt-smcopt.target); % Adjust the scale parameter
        scl(i) = smcopt.scale; % Scale parameter (for the jumping distribution in MH mutation step).
        mu = particles*weights; % Weighted average of the particles.
        z = particles-mu;
        R = z*(z'.*weights); % Weighted covariance matrix of the particles.
        t0 = tic;
        acpt_ = false(smcopt.particles, 1);
        tlogpostkernel = tlogpostkernel + (phi(i)-phi(i-1))*loglikelihood;
        [acpt_, particles, loglikelihood, tlogpostkernel] = ...
            randomwalk(funobj, chol(R, 'lower'), mu, scl(i), phi(i), acpt_, Prior, particles, loglikelihood, tlogpostkernel);
        smcopt.acpt = sum(acpt_)/smcopt.particles; % Acceptance rate.
        tt = toc(t0);
        acpt(i) = smcopt.acpt;
        if resampled_particle_swarm
            dprintf('%3u          %5.4f     %9.5E         %5.4f        %4.3f     %3s       %5.2f', i, phi(i), ESS(i), acpt(i), scl(i), 'yes', tt)
        else
            dprintf('%3u          %5.4f     %9.5E         %5.4f        %4.3f     %3s       %5.2f', i, phi(i), ESS(i), acpt(i), scl(i), 'no', tt)
        end
        if i==smcopt.steps
            iresample = kitagawa(weights);
            particles = particles(:,iresample);
        end
        save(sprintf('%s%sparticles-%u-%u.mat', SimulationFolder, filesep(), i, smcopt.steps), 'particles', 'tlogpostkernel', 'loglikelihood')
        resampled_particle_swarm = false;
    end

end

function [accept, particles, loglikelihood, tlogpostkernel] = randomwalk(funobj, RR, mu, scale, phi, accept, Prior, particles, loglikelihood, tlogpostkernel)

    parfor j=1:size(particles, 2)
        notvalid= true;
        while notvalid
            candidate = particles(:,j) + scale*(RR*randn(size(mu)));
            if Prior.admissible(candidate)
                [tlogpost, loglik] = tempered_likelihood(funobj, candidate, phi, Prior);
                if isfinite(loglik)
                    notvalid = false;
                    if rand<exp(tlogpost-tlogpostkernel(j))
                        accept(j) = true ;
                        particles(:,j) = candidate;
                        loglikelihood(j) = loglik;
                        tlogpostkernel(j) = tlogpost;
                    end
                end
            end
        end
    end
end
