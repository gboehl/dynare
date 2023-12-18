function display(pdraws, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, name, options_ident)
% display(pdraws, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, name, options_ident)
% -------------------------------------------------------------------------
% This function displays all identification analysis to the command line
% =========================================================================
% INPUTS
%   pdraws:             [SampleSize by totparam_nbr] parameter draws
%   ide_reducedform:    [structure] Containing results from identification
%                       analysis based on the reduced-form solution (Ratto
%                       and Iskrev, 2011).
%   ide_moments:        [structure] Containing results from identification
%                       analysis based on moments (Iskrev, 2010).
%   ide_spectrum:       [structure] Containing results from identification
%                       analysis based on the spectrum (Qu and Tkachenko, 2012).
%   ide_minimal:        [structure] Containing results from identification
%                       analysis based on the minimal state space system
%                       (Komunjer and Ng, 2011).
%   name:               [totparam_nbr by 1] string cell of parameter names
%   options_ident:      [structure] identification options
%   error_indicator     [structure] 
%                       boolean information on errors (1 is an error, 0 is no error)
%                       while computing the criteria

% -------------------------------------------------------------------------
% OUTPUTS
%   * all output is printed on the command line
% -------------------------------------------------------------------------
% This function is called by
%   * identification.run
% =========================================================================
% Copyright © 2010-2021 Dynare Team
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
[SampleSize, totparam_nbr] = size(pdraws);
no_identification_reducedform      = options_ident.no_identification_reducedform;
no_identification_moments          = options_ident.no_identification_moments;
no_identification_spectrum         = options_ident.no_identification_spectrum;
no_identification_minimal          = options_ident.no_identification_minimal;
tol_rank           = options_ident.tol_rank;
checks_via_subsets = options_ident.checks_via_subsets;

%% Display settings
disp('  '),
fprintf('Note that differences in the criteria could be due to numerical settings,\n')
fprintf('numerical errors or the method used to find problematic parameter sets.\n')
fprintf('Settings:\n')
if options_ident.analytic_derivation_mode == 0
    fprintf('    Derivation mode for Jacobians:                         Analytic using sylvester equations\n');
elseif options_ident.analytic_derivation_mode == 1
    fprintf('    Derivation mode for Jacobians:                         Analytic using kronecker products\n');
elseif options_ident.analytic_derivation_mode < 0
    fprintf('    Derivation mode for Jacobians:                         Numerical\n');
end
if checks_via_subsets
    fprintf('    Method to find problematic parameters:                 Rank condition on all possible subsets\n');
else
    fprintf('    Method to find problematic parameters:                 Nullspace and multicorrelation coefficients\n');
end
if options_ident.normalize_jacobians == 1
    fprintf('    Normalize Jacobians:                                   Yes\n');
else
    fprintf('    Normalize Jacobians:                                   No\n');
end
fprintf('    Tolerance level for rank computations:                 %s\n',num2str(options_ident.tol_rank));
fprintf('    Tolerance level for selecting nonzero columns:         %.0d\n',options_ident.tol_deriv);
fprintf('    Tolerance level for selecting nonzero singular values: %.0d\n',options_ident.tol_sv);


%% Display problematic parameter sets for different criteria in a loop
for jide = 1:4
    no_warning_message_display = 1;
    %% Set output strings depending on test
    if jide == 1
        strTest = 'REDUCED-FORM'; strJacobian = 'Tau'; strMeaning = 'Jacobian of steady state and reduced-form solution matrices';
        if ~no_identification_reducedform && ~ isempty(fieldnames(ide_reducedform))
            noidentification = 0; ide = ide_reducedform;
            if SampleSize == 1
                Jacob = ide.dREDUCEDFORM;
            end
        else %skip test
            noidentification = 1; no_warning_message_display = 0;
        end
    elseif jide == 2
        strTest = 'MINIMAL SYSTEM (Komunjer and Ng, 2011)'; strJacobian = 'Deltabar'; strMeaning = 'Jacobian of steady state and minimal system';
        if options_ident.order == 2
            strMeaning = 'Jacobian of first-order minimal system and second-order accurate mean';
        elseif options_ident.order == 3
            strMeaning = 'Jacobian of first-order minimal system and third-order accurate mean';
        end
        if ~no_identification_minimal && ~(length(ide_minimal.minimal_state_space)==1 && ide_minimal.minimal_state_space==0) && isfield(ide_minimal,'dMINIMAL')
            noidentification = 0; ide = ide_minimal;
            if SampleSize == 1
                Jacob = ide.dMINIMAL;
            end
        else %skip test
            noidentification = 1; no_warning_message_display = 0;
        end
    elseif jide == 3
        strTest = 'SPECTRUM (Qu and Tkachenko, 2012)'; strJacobian = 'Gbar'; strMeaning = 'Jacobian of mean and spectrum';
        if options_ident.order > 1
            strTest = 'SPECTRUM (Mutschler, 2015)';
        end
        if ~no_identification_spectrum && ~isempty(fieldnames(ide_spectrum))
            noidentification = 0; ide = ide_spectrum;
            if SampleSize == 1
                Jacob = ide.dSPECTRUM;
            end
        else %skip test
            noidentification = 1; no_warning_message_display = 0;
        end
    elseif jide == 4
        strTest = 'MOMENTS (Iskrev, 2010)'; strJacobian = 'J'; strMeaning = 'Jacobian of first two moments';
        if options_ident.order > 1
            strTest = 'MOMENTS (Mutschler, 2015)'; strJacobian = 'Mbar';
        end
        if ~no_identification_moments && ~isempty(fieldnames(ide_moments))
            noidentification = 0; ide = ide_moments;
            if SampleSize == 1
                Jacob = ide.si_dMOMENTS;
            end
        else %skip test
            noidentification = 1; no_warning_message_display = 0;
        end
    end

    if ~noidentification
        %% display problematic parameters computed by identifcation_checks.m
        if ~checks_via_subsets
            EffectiveSampleSize=SampleSize;
            non_minimal_state_space_error=0;
            if SampleSize>1 && jide==2  && any(~ide.minimal_state_space)
                EffectiveSampleSize=SampleSize-sum(~ide.minimal_state_space);
                non_minimal_state_space_error=1;
                ide.ino(~ide.minimal_state_space,:)=[];
                ide.ind0(~ide.minimal_state_space,:)=[];
                ide.jweak_pair(~ide.minimal_state_space,:)=[];
                ide.minimal_state_space(~ide.minimal_state_space,:)=[];
            end
            if any(ide.ino) || any(any(ide.ind0==0)) || any(any(ide.jweak_pair)) || non_minimal_state_space_error
                no_warning_message_display=0;
                skipline()
                disp([upper(strTest), ':'])
                disp('    !!!WARNING!!!');
                if SampleSize>1
                    if non_minimal_state_space_error
                        fprintf('\n    The minimal state space could not be computed for %u out of %u cases.\n',SampleSize-EffectiveSampleSize,SampleSize);
                    end
                    if jide==2
                        if sum(ide.ino & ide.minimal_state_space)>0
                        disp(['    The rank of ', strJacobian, ' (', strMeaning, ') is deficient for ', num2str(sum(ide.ino & ide.minimal_state_space)),' out of ',int2str(EffectiveSampleSize),' effective MC runs!'  ])
                        end
                    else
                        disp(['    The rank of ', strJacobian, ' (', strMeaning, ') is deficient for ', num2str(sum(ide.ino~=0)),' out of ',int2str(EffectiveSampleSize),' effective MC runs!'  ]),
                    end
                else
                    disp(['    The rank of ', strJacobian, ' (', strMeaning, ') is deficient!']),
                end
                skipline()
                for j=1:totparam_nbr
                    if any(ide.ind0(:,j)==0)
                        pno = 100*length(find(ide.ind0(:,j)==0))/EffectiveSampleSize;
                        if SampleSize>1
                            disp(['    ',name{j},' is not identified for ',num2str(pno),'% of MC runs!' ])
                        else
                            disp(['    ',name{j},' is not identified!' ])
                        end
                    end
                end
                npairs=size(ide.jweak_pair,2);
                jmap_pair=dyn_unvech(1:npairs);
                jstore=[];
                for j=1:npairs
                    iweak = length(find(ide.jweak_pair(:,j)));
                    if iweak
                        [jx,jy]=find(jmap_pair==j);
                        jstore=[jstore jx(1) jy(1)];
                        if SampleSize > 1
                            disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear for ',num2str((iweak)/EffectiveSampleSize*100),'% of MC runs!' ])
                        else
                            disp(['    [',name{jx(1)},',',name{jy(1)},'] are PAIRWISE collinear!' ])
                        end
                    end
                end
                for j=1:totparam_nbr
                    iweak = length(find(ide.jweak(:,j)));
                    if iweak && ~ismember(j,jstore)
                        if SampleSize>1
                            disp(['    ',name{j},' is collinear w.r.t. all other parameters for ',num2str(iweak/EffectiveSampleSize*100),'% of MC runs!' ])
                        else
                            disp(['    ',name{j},' is collinear w.r.t. all other parameters!' ])
                        end
                    end
                end
            end

            %% display problematic parameters computed by identification.checks_via_subsets
        elseif checks_via_subsets
            if ide.rank < size(Jacob,2)
                no_warning_message_display = 0;
                skipline()
                disp([upper(strTest), ':'])
                disp('    !!!WARNING!!!');
                if SampleSize>1
                    disp(['    The rank of ', strJacobian, ' (', strMeaning, ') is deficient for ', num2str(length(find(ide.ino))),' out of ',int2str(SampleSize),' MC runs!'  ]),
                else
                    disp(['    The rank of ', strJacobian, ' (', strMeaning, ') is deficient!']),
                end
                if all(cellfun(@isempty,ide.problpars))
                    disp(['    No problematic parameter combinations with maximum dimension ', num2str(size(ide.problpars,2)), ' were found. Increase max_dim_subsets_groups.']);
                    skipline()
                else
                    disp(['    Displaying problematic parameter combinations (with maximum dimension ', num2str(size(ide.problpars,2)), '):']);
                    skipline()
                    probparamset_nbr = 0;
                    for jset = 1:size(ide.problpars,2)
                        if isempty(ide.problpars{jset}) == 0
                            for jrow = 1:size(ide.problpars{jset},1)
                                for j = transpose(ide.problpars{jset}(jrow,:))
                                    probparamset_nbr = probparamset_nbr + 1;
                                    %pno = 100*length(find(ide.ind0(:,j)==0))/SampleSize;
                                    problparnamestring = strjoin(eval(['[', sprintf('name(%d), ', j), ']']),',');
                                    if SampleSize > 1
                                        if length(j) == 1
                                            disp(['    ',problparnamestring,' is not identified for ',num2str(pno),'% of MC runs!' ])
                                        else
                                            disp(['    [',problparnamestring,'] are collinear (with tol = ', num2str(tol_rank), ') for ',num2str((iweak)/SampleSize*100),'% of MC runs!' ])
                                        end
                                    else
                                        if length(j) == 1
                                            disp(['    ',problparnamestring, ' is not identified!' ])
                                        else
                                            disp(['    [',problparnamestring, '] are collinear!' ])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %% All parameters are identified
    if no_warning_message_display
        skipline()
        disp([upper(strTest), ':']);
        disp(['    All parameters are identified in the ', strMeaning, ' (rank(', strJacobian, ') is full with tol = ', num2str(tol_rank), ').' ]),
    end
end



%% Advanced identification patterns
if SampleSize==1 && options_ident.advanced && ~isempty(fieldnames(ide_moments))
    skipline()
    for j=1:size(ide_moments.cosndMOMENTS,2)
        pax=NaN(totparam_nbr,totparam_nbr);
        fprintf('\n')
        disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
        fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
        for i=1:totparam_nbr
            namx='';
            for in=1:j
                dumpindx = ide_moments.pars{i,j}(in);
                if isnan(dumpindx)
                    namx=[namx ' ' sprintf('%-15s','--')];
                else
                    namx=[namx ' ' sprintf('%-15s',name{dumpindx})];
                    pax(i,dumpindx)=ide_moments.cosndMOMENTS(i,j);
                end
            end
            fprintf('%-15s [%s] %14.7f\n',name{i},namx,ide_moments.cosndMOMENTS(i,j))
        end
    end
end
