function [nvar,vartan,CovarFileNumber] = dsge_simulated_theoretical_covariance(SampleSize,nar,M_,options_,oo_,type)
% function [nvar,vartan,CovarFileNumber] = dsge_simulated_theoretical_covariance(SampleSize,nar,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the endogenous
% variables second order moments.
%
% INPUTS
%   SampleSize   [integer]       scalar, number of simulations.
%   M_           [structure]     Dynare structure describing the model.
%   options_     [structure]     Dynare structure defining global options.
%   oo_          [structure]     Dynare structure where the results are saved.
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   nvar              [integer]  nvar is the number of stationary variables.
%   vartan            [char]     array of characters (with nvar rows).
%   CovarFileNumber   [integer]  scalar, number of prior or posterior data files (for covariance).

% Copyright (C) 2007-2021 Dynare Team
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

nodecomposition = 1;

% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    NumberOfDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]));
    posterior = 1;
elseif strcmpi(type,'prior')
    NumberOfDrawsFiles = length(dir([M_.dname '/prior/draws/' type '_draws*' ]));
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    disp('dsge_simulated_theoretical_covariance:: Unknown type!')
    error();
end

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments*'])
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations*']);
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_Prior2ndOrderMoments*'])
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PriorCorrelations*']);
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
endo_names=options_.varlist;
[ivar,vartan] = get_variables_list(options_,M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

if options_.pruning
    obs_var=NaN(nvar,1);
    for i=1:nvar
        obs_var(i,1) = find(strcmp(M_.endo_names(ivar(i),:), M_.endo_names(oo_.dr.order_var)));
    end
end

options_.ar = nar;

% Number of lines in posterior data files.
MaXNumberOfCovarLines = ceil(options_.MaxNumberOfBytes/(nvar*(nvar+1)/2)/8);
MaXNumberOfCorrLines = ceil(options_.MaxNumberOfBytes/(nvar*nvar*nar)/8);

if SampleSize<=MaXNumberOfCovarLines
    Covariance_matrix = zeros(SampleSize,nvar*(nvar+1)/2);
    NumberOfCovarFiles = 1;
else
    Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
    NumberOfLinesInTheLastCovarFile = mod(SampleSize,MaXNumberOfCovarLines);
    NumberOfCovarFiles = ceil(SampleSize/MaXNumberOfCovarLines);
end
if SampleSize<=MaXNumberOfCorrLines
    Correlation_array = zeros(SampleSize,nvar,nvar,nar);
    NumberOfCorrFiles = 1;
else
    Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
    NumberOfLinesInTheLastCorrFile = mod(SampleSize,MaXNumberOfCorrLines);
    NumberOfCorrFiles = ceil(SampleSize/MaXNumberOfCorrLines);
end

NumberOfCovarLines = rows(Covariance_matrix);
CovarFileNumber = 1;
NumberOfCorrLines = rows(Correlation_array);
CorrFileNumber = 1;

% Compute 2nd order moments and save them in *_[Posterior, Prior]2ndOrderMoments* files
linea_cov = 0;
linea_corr = 0;
for file = 1:NumberOfDrawsFiles
    if posterior
        temp=load([M_.dname '/metropolis/' M_.fname '_' type '_draws' num2str(file) ]);
    else
        temp=load([M_.dname '/prior/draws/' type '_draws' num2str(file) ]);
    end
    NumberOfDraws = rows(temp.pdraws);
    isdrsaved = columns(temp.pdraws)-1;
    for linee = 1:NumberOfDraws
        linea_cov = linea_cov+1;
        linea_corr = linea_corr+1;
        if isdrsaved
            M_=set_parameters_locally(M_,temp.pdraws{linee,1});% Needed to update the covariance matrix of the state innovations.
            dr = temp.pdraws{linee,2};
        else
            M_=set_parameters_locally(M_,temp.pdraws{linee,1});
            [dr,info,M_,oo_] = compute_decision_rules(M_,options_,oo_);
        end
        if ~options_.pruning
            tmp = th_autocovariances(dr,ivar,M_,options_,nodecomposition);
        else
            pruned_state_space = pruned_state_space_system(M_, options_, dr, obs_var, options_.ar, 1, 0);
            tmp{1} = pruned_state_space.Var_y;            
            for i=1:nar
                tmp{i+1} = pruned_state_space.Corr_yi(:,:,i);                
            end
        end
        for i=1:nvar
            for j=i:nvar
                Covariance_matrix(linea_cov,symmetric_matrix_index(i,j,nvar)) = tmp{1}(i,j);
            end
        end
        for i=1:nar
            Correlation_array(linea_corr,:,:,i) = tmp{i+1};
        end

        if linea_cov == NumberOfCovarLines
            if posterior
                save([ M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments' int2str(CovarFileNumber) '.mat' ],'Covariance_matrix','endo_names');
            else
                save([ M_.dname '/prior/moments/' M_.fname '_Prior2ndOrderMoments' int2str(CovarFileNumber) '.mat' ],'Covariance_matrix','endo_names');
            end
            CovarFileNumber = CovarFileNumber + 1;
            linea_cov = 0;
            test = CovarFileNumber-NumberOfCovarFiles;
            if ~test% Prepare the last round...
                Covariance_matrix = zeros(NumberOfLinesInTheLastCovarFile,nvar*(nvar+1)/2);
                NumberOfCovarLines = NumberOfLinesInTheLastCovarFile;
            elseif test<0
                Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
            else
                clear('Covariance_matrix');
            end
        end
        if linea_corr == NumberOfCorrLines
            if posterior
                save([ M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations' int2str(CorrFileNumber) '.mat' ],'Correlation_array','endo_names');
            else
                save([ M_.dname '/prior/moments/' M_.fname '_PriorCorrelations' int2str(CorrFileNumber) '.mat' ],'Correlation_array','endo_names');
            end
            CorrFileNumber = CorrFileNumber + 1;
            linea_corr = 0;
            test = CorrFileNumber-NumberOfCorrFiles;
            if ~test% Prepare the last round...
                Correlation_array = zeros(NumberOfLinesInTheLastCorrFile,nvar,nvar,nar);
                NumberOfCorrLines = NumberOfLinesInTheLastCorrFile;
                CorrFileNumber = CorrFileNumber - 1;
            elseif test<0
                Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
            else
                clear('Correlation_array');
            end
        end        
    end
end
