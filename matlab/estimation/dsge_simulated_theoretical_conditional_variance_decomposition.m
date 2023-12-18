function [nvar,vartan,NumberOfConditionalDecompFiles] = ...
    dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,Steps,M_,options_,oo_,type)
% function [nvar,vartan,NumberOfConditionalDecompFiles] = ...
%     dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,Steps,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the conditional variance
% decomposition of the endogenous variables (or a subset of the endogenous variables).
%
% INPUTS
%   SampleSize   [integer]       scalar, number of simulations.
%   Steps        [integers]      horizons at which to conduct decomposition
%   M_           [structure]     Dynare structure describing the model.
%   options_     [structure]     Dynare structure defining global options.
%   oo_          [structure]     Dynare structure where the results are saved.
%   type         [string]        'prior' or 'posterior'
%
%
% OUTPUTS
%   nvar                             [integer]  nvar is the number of stationary variables.
%   vartan                           [char]     array of characters (with nvar rows).
%   NumberOfConditionalDecompFiles   [integer]  scalar, number of prior or posterior data files (for covariance).

% Copyright © 2009-2023 Dynare Team
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


% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    NumberOfDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]));
    posterior = 1;
elseif strcmpi(type,'prior')
    NumberOfDrawsFiles = length(dir([M_.dname '/prior/draws/' type '_draws*' ]));
    CheckPath('prior/moments',M_.dname);
    posterior = 0;
else
    error('dsge_simulated_theoretical_conditional_variance_decomposition:: Unknown type requested!')
end

%delete old stale files before creating new ones
if posterior
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecomposition*'])
    delete_stale_file([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecompositionME*'])
else
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecomposition*'])
    delete_stale_file([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecompositionME*'])
end

% Set varlist (vartan)
if ~posterior
    if isfield(options_,'varlist')
        temp = options_.varlist;
    end
    options_.varlist = options_.prior_analysis_endo_var_list;
end
endo_names=options_.varlist;
[ivar,vartan ] = get_variables_list(options_, M_);
if ~posterior
    if exist('temp','var')
        options_.varlist = temp;
    end
end
nvar = length(ivar);

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;

NumberOfSavedElementsPerSimulation = nvar*M_.exo_nbr*length(Steps);
MaXNumberOfConditionalDecompLines = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation/8);

[ME_present,observable_pos_requested_vars] = check_measurement_error_requested_vars(M_,options_,ivar);

if ME_present && ~isempty(observable_pos_requested_vars)
    nobs_ME=length(observable_pos_requested_vars);
    NumberOfSavedElementsPerSimulation_ME = nobs_ME*(M_.exo_nbr+1)*length(Steps);
    MaXNumberOfConditionalDecompLines_ME = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation_ME/8);
end

if SampleSize<=MaXNumberOfConditionalDecompLines
    Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,SampleSize);
    NumberOfConditionalDecompFiles = 1;
else
    Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,MaXNumberOfConditionalDecompLines);
    NumberOfLinesInTheLastConditionalDecompFile = mod(SampleSize,MaXNumberOfConditionalDecompLines);
    NumberOfConditionalDecompFiles = ceil(SampleSize/MaXNumberOfConditionalDecompLines);
end

if ME_present
    if SampleSize<=MaXNumberOfConditionalDecompLines_ME
        Conditional_decomposition_array_ME = zeros(nobs_ME,length(Steps),M_.exo_nbr+1,SampleSize);
        NumberOfConditionalDecompFiles_ME = 1;
    else
        Conditional_decomposition_array_ME = zeros(nobs_ME,length(Steps),M_.exo_nbr+1,SampleSize);
        NumberOfLinesInTheLastConditionalDecompFile_ME = mod(SampleSize,MaXNumberOfConditionalDecompLines_ME);
        NumberOfConditionalDecompFiles_ME = ceil(SampleSize/MaXNumberOfConditionalDecompLines_ME);
    end
    NumberOfConditionalDecompLines_ME = size(Conditional_decomposition_array_ME,4);
    ConditionalDecompFileNumber_ME = 0;
end

NumberOfConditionalDecompLines = size(Conditional_decomposition_array,4);
ConditionalDecompFileNumber = 0;

linea = 0;
linea_ME = 0;
for file = 1:NumberOfDrawsFiles
    if posterior
        temp=load([M_.dname '/metropolis/' M_.fname '_' type '_draws' num2str(file) ]);
    else
        temp=load([M_.dname '/prior/draws/' type '_draws' num2str(file) ]);
    end
    isdrsaved = columns(temp.pdraws)-1;
    NumberOfDraws = rows(temp.pdraws);
    for linee = 1:NumberOfDraws
        linea = linea+1;
        linea_ME = linea_ME+1;
        if isdrsaved
            M_=set_parameters_locally(M_,temp.pdraws{linee,1});% Needed to update the covariance matrix of the state innovations.
            dr = temp.pdraws{linee,2};
        else
            M_=set_parameters_locally(M_,temp.pdraws{linee,1});
			[dr,~,M_.params] = compute_decision_rules(M_,options_,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        end

        M_ = set_measurement_errors(temp.pdraws{linee,1},temp.estim_params_,M_);
        [ConditionalVarianceDecomposition, ConditionalVarianceDecomposition_ME] = ...
        conditional_variance_decomposition(M_,options_,dr, Steps, ivar);
        Conditional_decomposition_array(:,:,:,linea) =ConditionalVarianceDecomposition;
        if ME_present
            Conditional_decomposition_array_ME(:,:,:,linea) =ConditionalVarianceDecomposition_ME;
        end
        if linea == NumberOfConditionalDecompLines
            ConditionalDecompFileNumber = ConditionalDecompFileNumber + 1;
            linea = 0;
            if posterior
                save([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecomposition' int2str(ConditionalDecompFileNumber) '.mat' ], ...
                     'Conditional_decomposition_array','endo_names');
            else
                save([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecomposition' int2str(ConditionalDecompFileNumber) '.mat' ], ...
                     'Conditional_decomposition_array','endo_names');
            end
            if (ConditionalDecompFileNumber==NumberOfConditionalDecompFiles-1)% Prepare last round.
                Conditional_decomposition_array = zeros(nvar, length(Steps),M_.exo_nbr,NumberOfLinesInTheLastConditionalDecompFile) ;
                NumberOfConditionalDecompLines = NumberOfLinesInTheLastConditionalDecompFile;
            elseif ConditionalDecompFileNumber<NumberOfConditionalDecompFiles-1
                Conditional_decomposition_array = zeros(nvar,length(Steps),M_.exo_nbr,MaXNumberOfConditionalDecompLines);
            else
                clear('Conditional_decomposition_array');
            end
        end
        %with measurement error
        if ME_present
            if linea_ME == NumberOfConditionalDecompLines_ME
                ConditionalDecompFileNumber_ME = ConditionalDecompFileNumber_ME + 1;
                linea_ME = 0;
                if posterior
                    save([M_.dname '/metropolis/' M_.fname '_PosteriorConditionalVarianceDecompME' int2str(ConditionalDecompFileNumber_ME) '.mat' ], ...
                         'Conditional_decomposition_array_ME','endo_names');
                else
                    save([M_.dname '/prior/moments/' M_.fname '_PriorConditionalVarianceDecompME' int2str(ConditionalDecompFileNumber_ME) '.mat' ], ...
                         'Conditional_decomposition_array_ME','endo_names');
                end
                if (ConditionalDecompFileNumber_ME==NumberOfConditionalDecompFiles_ME-1)% Prepare last round.
                    Conditional_decomposition_array_ME = zeros(nobs_ME, length(Steps),M_.exo_nbr+1,NumberOfLinesInTheLastConditionalDecompFile_ME) ;
                    NumberOfConditionalDecompLines_ME = NumberOfLinesInTheLastConditionalDecompFile_ME;
                elseif ConditionalDecompFileNumber_ME<NumberOfConditionalDecompFiles_ME-1
                    Conditional_decomposition_array_ME = zeros(nobs_ME,length(Steps),M_.exo_nbr+1,MaXNumberOfConditionalDecompLines_ME);
                else
                    clear('Conditional_decomposition_array_ME');
                end
            end
        end
    end
end

options_.ar = nar;
