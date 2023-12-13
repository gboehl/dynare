function [mean, covariance, mode, kernel_at_the_mode] = compute_mh_covariance_matrix(names, fname, dname, outputFolderName)

% Estimation of the posterior covariance matrix, posterior mean, posterior mode and evaluation of the posterior kernel at the
% estimated mode, using posterior draws from a metropolis-hastings.
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

% Copyright © 2006-2023 Dynare Team
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

MetropolisFolder = CheckPath('metropolis',dname);
BaseName = [MetropolisFolder filesep fname];

record=load_last_mh_history_file(MetropolisFolder, fname);

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine   = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));

[nblck, n] = size(record.LastParameters);

kernel_at_the_mode = -Inf;
mean = zeros(n,1);
mode = NaN(n,1);
covariance = zeros(n,n);
offset = 0;

for b=1:nblck
    first_line = FirstLine;
    for n = FirstMhFile:TotalNumberOfMhFiles
        load([ BaseName '_mh' int2str(n) '_blck' int2str(b) '.mat'],'x2','logpo2');
        [tmp, idx] = max(logpo2);
        if tmp>kernel_at_the_mode
            kernel_at_the_mode = tmp;
            mode = x2(idx,:);
        end
        [mean, covariance, offset] = recursive_moments(mean, covariance, x2(first_line:end,:), offset);
        first_line = 1;
    end
end

mode = transpose(mode);
