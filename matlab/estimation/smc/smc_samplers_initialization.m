function [particles, tlogpostkernel, loglikelihood, SimulationFolder] = smc_samplers_initialization(funobj, sampler, n, Prior, SimulationFolder, nsteps)

% Initialize SMC samplers by drawing initial particles in the prior distribution.
%
% INPUTS
% - TargetFun        [char]     string specifying the name of the objective function (posterior kernel).
% - sampler          [char]     name of the sampler.
% - n                [integer]  scalar, number of particles.
% - mh_bounds        [double]   p×2 matrix defining lower and upper bounds for the estimated parameters.
% - dataset_         [dseries]  sample
% - dataset_info     [struct]   informations about the dataset
% - options_         [struct]   dynare's options
% - M_               [struct]   model description
% - estim_params_    [struct]   estimated parameters
% - bayestopt_       [struct]   estimated parameters
% - oo_              [struct]   outputs
%
% OUTPUTS
% - ix2                   [double]    p×n matrix of particles
% - ilogpo2               [double]    n×1 vector of posterior kernel values for the particles
% - iloglik2              [double]    n×1 vector of likelihood values for the particles
% - ModelName             [string]    name of the mod-file
% - MetropolisFolder      [string]    path to the Metropolis subfolder
% - bayestopt_            [structure] estimation options structure
%
% SPECIAL REQUIREMENTS
%   None.

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

dprintf('Estimation:%s: Initialization...', sampler)

% Delete old mat files storign particles if any...
matfiles = sprintf('%s%sparticles*.mat', SimulationFolder, filesep());
files = dir(matfiles);
if ~isempty(files)
    delete(matfiles);
    dprintf('Estimation:%s: Old %s-files successfully erased.', sampler, sampler)
end

% Simulate a pool of particles characterizing the prior distribution (with the additional constraint that the likelihood is finite)
set_dynare_seed('default');
dprintf('Estimation:%s: Searching for initial values...', sampler);
particles = zeros(Prior.length(), n);
tlogpostkernel = zeros(n, 1);
loglikelihood = zeros(n, 1);

t0 = tic;
parfor j=1:n
    notvalid = true;
    while notvalid
        candidate = Prior.draw();
        if Prior.admissible(candidate)
            particles(:,j) = candidate;
            [tlogpostkernel(j), loglikelihood(j)] = tempered_likelihood(funobj, candidate, 0.0, Prior);
            if isfinite(loglikelihood(j)) % if returned log-density is Inf or Nan (penalized value)
                notvalid = false;
            end
        end
    end
end
tt = toc(t0);

save(sprintf('%s%sparticles-1-%u.mat', SimulationFolder, filesep(), nsteps), 'particles', 'tlogpostkernel', 'loglikelihood')
dprintf('Estimation:%s: Initial values found (%.2fs)', sampler, tt)
skipline()
