function [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_, bayestopt_, outputFolderName)
% function [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_, bayestopt_, outputFolderName)
% Computes the marginal density
%
% INPUTS
%   options_         [structure]    Dynare options structure
%   estim_params_    [structure]    Dynare estimation parameter structure
%   M_               [structure]    Dynare model structure
%   oo_              [structure]    Dynare results structure
%   outputFolderName [string]       name of folder with results
%
% OUTPUTS
%   marginal:        [double]       marginal density (modified harmonic mean)
%   oo_              [structure]    Dynare results structure
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2005-2023 Dynare Team
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
if nargin < 6
    outputFolderName = 'Output';
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

record=load_last_mh_history_file(MetropolisFolder, ModelName);

[nblck, npar] = size(record.LastParameters);

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
TODROP = floor(options_.mh_drop*TotalNumberOfMhDraws);

fprintf('marginal density: I''m computing the posterior mean and covariance... ');
[posterior_mean, posterior_covariance, posterior_mode, posterior_kernel_at_the_mode] = compute_posterior_covariance_matrix(bayestopt_.name, M_.fname, M_.dname, options_, outputFolderName);

MU = transpose(posterior_mean);
SIGMA = posterior_covariance;
lpost_mode = posterior_kernel_at_the_mode;
xparam1 = posterior_mean;
hh = inv(SIGMA);
fprintf(' Done!\n');
if ~isfield(oo_,'posterior_mode') || (options_.mh_replic && isequal(options_.posterior_sampler_options.posterior_sampling_method,'slice'))
    oo_=fill_mh_mode(posterior_mode',NaN(npar,1),M_,options_,estim_params_,oo_,'posterior');
end

% save the posterior mean and the inverse of the covariance matrix
% (usefull if the user wants to perform some computations using
% the posterior mean instead of the posterior mode ==> ).
parameter_names = bayestopt_.name;
save([M_.dname filesep outputFolderName filesep M_.fname '_mean.mat'],'xparam1','hh','parameter_names','SIGMA');

fprintf('marginal density: I''m computing the posterior log marginal density (modified harmonic mean)... ');
try 
    % use this robust option to avoid inf/nan
    logdetSIGMA = 2*sum(log(diag(chol(SIGMA)))); 
catch
    % in case SIGMA is not positive definite
    logdetSIGMA = nan;
    fprintf('marginal density: the covariance of MCMC draws is not positive definite. You may have too few MCMC draws.');
end
invSIGMA = hh;
marginal = zeros(9,2);
linee = 0;
check_coverage = 1;
increase = 1;
while check_coverage
    for p = 0.1:0.1:0.9
        critval = chi2inv(p,npar);
        ifil = FirstLine;
        tmp = 0;
        for n = FirstMhFile:TotalNumberOfMhFiles
            for b=1:nblck
                load([ BaseName '_mh' int2str(n) '_blck' int2str(b) '.mat'],'x2','logpo2');
                EndOfFile = size(x2,1);
                for i = ifil:EndOfFile
                    deviation  = ((x2(i,:)-MU)*invSIGMA*(x2(i,:)-MU)')/increase;
                    if deviation <= critval
                        lftheta = -log(p)-(npar*log(2*pi)+(npar*log(increase)+logdetSIGMA)+deviation)/2;
                        tmp = tmp + exp(lftheta - logpo2(i) + lpost_mode);
                    end
                end
            end
            ifil = 1;
        end
        linee = linee + 1;
        warning_old_state = warning;
        warning off;
        marginal(linee,:) = [p, lpost_mode-log(tmp/((TotalNumberOfMhDraws-TODROP)*nblck))];
        warning(warning_old_state);
    end
    if abs((marginal(9,2)-marginal(1,2))/marginal(9,2)) > options_.marginal_data_density.harmonic_mean.tolerance || isinf(marginal(1,2))
        fprintf('\n')
        if increase == 1
            disp('marginal density: The support of the weighting density function is not large enough...')
            disp('marginal density: I increase the variance of this distribution.')
            increase = 1.2*increase;
            linee    = 0;
        else
            disp('marginal density: Let me try again.')
            increase = 1.2*increase;
            linee    = 0;
            if increase > 20
                check_coverage = 0;
                clear invSIGMA detSIGMA increase;
                disp('marginal density: There''s probably a problem with the modified harmonic mean estimator.')
            end
        end
    else
        check_coverage = 0;
        clear invSIGMA detSIGMA increase;
        fprintf('Done!\n')
    end
end

oo_.MarginalDensity.ModifiedHarmonicMean = mean(marginal(:,2));

