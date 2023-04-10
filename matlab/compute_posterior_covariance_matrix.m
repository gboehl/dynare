function [mu, covariance, mode, kernel_at_the_mode] = compute_posterior_covariance_matrix(names, fname, dname, options_, outputFolderName)

% Estimation of the posterior covariance matrix, posterior mean, posterior mode and evaluation of the posterior kernel at the
% estimated mode, using posterior draws from a metropolis-hastings. The estimated posterior mode and covariance matrix are saved in
% a file <fname>_mh_mode.mat, hssmc_mode.mat or dsmh__mode.mat under <dname>/<outputFolderName>/.
%
% INPUTS
% - names                    [cell]     n×1 cell array of row char arrays, names of the estimated parameters.
% - fname                    [char]     name of the model
% - dname                    [char]     name of subfolder with output files
% - outputFolderName         [char]     name of directory to store results
%
% OUTPUTS
% - mean                     [double]   n×1 vector, posterior expectation of the parameters.
% - covariance               [double]   n×n matrix, posterior covariance of the parameters.
% - mode                     [double]   n×1 vector, posterior mode of the parameters.
% - kernel_at_the_mode       [double]   scalar, value of the posterior kernel at the mode.

% Copyright © 2023 Dynare Team
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


if nargin<5
    outputFolderName = 'Output';
end

if ishssmc(options_)
    % Load draws from the posterior distribution
    pfiles = dir(sprintf('%s/hssmc/particles-*.mat', dname));
    posterior = load(sprintf('%s/hssmc/particles-%u-%u.mat', dname, length(pfiles), length(pfiles)));
    % Get the posterior mode
    [kernel_at_the_mode, id] = max(posterior.tlogpostkernel);
    mode = posterior.particles(:,id);
    % Compute the posterior mean
    mu = sum(posterior.particles, 2)/length(posterior.tlogpostkernel);
    % Compute the posterior covariance
    covariance = (posterior.particles-mu)*(posterior.particles-mu)'/length(posterior.tlogpostkernel);
else
    [mu, covariance, mode, kernel_at_the_mode] = compute_mh_covariance_matrix(names, fname, dname, outputFolderName);
end

xparam1 = mode;
hh = inv(covariance);
fval = kernel_at_the_mode;
parameter_names = names;

if ishssmc(options_)
    save(sprintf('%s/%s/hssmc_mode.mat', dname, outputFolderName), 'xparam1', 'hh', 'fval', 'parameter_names');
else
    save(sprintf('%s/%s/%s_mh_mode.mat', dname, outputFolderName, fname), 'xparam1', 'hh', 'fval', 'parameter_names');
end
