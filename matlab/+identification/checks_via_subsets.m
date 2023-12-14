function [ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = checks_via_subsets(ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident,error_indicator)
%[ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = checks_via_subsets(ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident,error_indicator)
% -------------------------------------------------------------------------
% Finds problematic sets of paramters via checking the necessary rank condition
% of the Jacobians for all possible combinations of parameters. The rank is
% computed via an inbuild function based on the SVD, similar to matlab's
% rank. The idea is that once we have the Jacobian for all parameters, we
% can easily set up Jacobians containing all combinations of parameters by
% picking the relevant columns/elements of the full Jacobian. Then the rank
% of these smaller Jacobians indicates whether this paramter combination is
% identified or not. To speed up computations:
% (1) single parameters are removed from possible higher-order sets,
% (2) for parameters that are collinear, i.e. rank failure for 2 element sets,
% we replace the second parameter by the first one, and then compute
% higher-order combinations [uncommented]
% (3) all lower-order problematic sets are removed from higher-order sets
% by an inbuild function
% (4) we could replace nchoosek by a mex version, e.g. VChooseK
% (https://de.mathworks.com/matlabcentral/fileexchange/26190-vchoosek) as
% nchoosek could be the bottleneck in terms of speed (and memory for models
% with totparam_nbr > 150)
% =========================================================================
% INPUTS
%   ide_reducedform:    [structure] Containing results from identification
%                       analysis based on the reduced-form solution (Ratto
%                       and Iskrev, 2011). If either options_ident.no_identification_reducedform
%                       or error_indicator.identification_reducedform are 1 then the search for problematic parameter sets will be skipped
%   ide_moments:        [structure] Containing results from identification
%                       analysis based on moments (Iskrev, 2010). If either
%                       options_ident.no_identification_moments or error_indicator.identification_moments are 1 then the search for
%                       problematic parameter sets will be skipped
%   ide_spectrum:       [structure] Containing results from identification
%                       analysis based on the spectrum (Qu and Tkachenko, 2012).
%                       If either options_ident.no_identification_spectrum or error_indicator.identification_spectrum are 1 then the search for
%                       problematic parameter sets will be skipped
%   ide_minimal:        [structure] Containing results from identification
%                       analysis based on the minimal state space system
%                       (Komunjer and Ng, 2011). If either options_ident.no_identification_minimal or
%                       error_indicator.identification_minimal are 1 then the search for problematic parameter sets will be skipped
%   totparam_nbr:       [integer] number of estimated stderr, corr and model parameters
%   numzerotolrank:     [double] tolerance level for rank compuations
%   error_indicator     [structure] indicators whether objects could be computed
% -------------------------------------------------------------------------
% OUTPUTS
%   ide_reducedform, ide_moments, ide_spectrum, ide_minimal are augmented by the
%   following fields:
%   * problpars:  [1 by totparam_nbr] cell with the following structure for j=1:totparam_nbr
%                  problpars{j}:   [nonidentified_j_set_parameters_nbr by j]
%                  matrix with j collinear parameters in each row
%   * rank:       [integer] rank of Jacobian
% -------------------------------------------------------------------------
% This function is called by
%   * identification.analysis.m
% =========================================================================
% Copyright © 2019-2021 Dynare Team
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

%% initialize output objects and get options
no_identification_dynamic             = 0; %always compute dynamic
no_identification_reducedform         = options_ident.no_identification_reducedform;
no_identification_moments             = options_ident.no_identification_moments;
no_identification_spectrum            = options_ident.no_identification_spectrum;
no_identification_minimal             = options_ident.no_identification_minimal;
tol_rank               = options_ident.tol_rank;
max_dim_subsets_groups = options_ident.max_dim_subsets_groups;
dynamic_problpars      = cell(1,max_dim_subsets_groups);
reducedform_problpars  = cell(1,max_dim_subsets_groups);
moments_problpars      = cell(1,max_dim_subsets_groups);
spectrum_problpars     = cell(1,max_dim_subsets_groups);
minimal_problpars      = cell(1,max_dim_subsets_groups);

indtotparam            = 1:totparam_nbr; %initialize index of parameters

%% Prepare Jacobians and check rank
% initialize linear rational expectations model
if ~no_identification_dynamic
    dDYNAMIC = ide_dynamic.dDYNAMIC;
    dDYNAMIC(ide_dynamic.ind_dDYNAMIC,:) = dDYNAMIC(ide_dynamic.ind_dDYNAMIC,:)./ide_dynamic.norm_dDYNAMIC; %normalize
    if totparam_nbr > modparam_nbr
        dDYNAMIC = [zeros(size(ide_dynamic.dDYNAMIC,1),totparam_nbr-modparam_nbr) dDYNAMIC]; %add derivatives wrt stderr and corr parameters
    end
    if strcmp(tol_rank,'robust')
        rank_dDYNAMIC = rank(full(dDYNAMIC)); %compute rank with imposed tolerance level
    else
        rank_dDYNAMIC = rank(full(dDYNAMIC),tol_rank); %compute rank with imposed tolerance level
    end
    ide_dynamic.rank = rank_dDYNAMIC;
    % check rank criteria for full Jacobian
    if rank_dDYNAMIC == totparam_nbr
        % all parameters are identifiable
        no_identification_dynamic = 1; %skip in the following
        indparam_dDYNAMIC = [];
    else
        % there is lack of identification
        indparam_dDYNAMIC = indtotparam; %initialize for nchoosek
    end
else
    indparam_dDYNAMIC = []; %empty for nchoosek
end

% initialize for reduced form solution criteria
if ~no_identification_reducedform && ~error_indicator.identification_reducedform
    dREDUCEDFORM = ide_reducedform.dREDUCEDFORM;
    dREDUCEDFORM(ide_reducedform.ind_dREDUCEDFORM,:) = dREDUCEDFORM(ide_reducedform.ind_dREDUCEDFORM,:)./ide_reducedform.norm_dREDUCEDFORM; %normalize
    if strcmp(tol_rank,'robust')
        rank_dREDUCEDFORM = rank(full(dREDUCEDFORM)); %compute rank with imposed tolerance level
    else
        rank_dREDUCEDFORM = rank(full(dREDUCEDFORM),tol_rank); %compute rank with imposed tolerance level
    end
    ide_reducedform.rank = rank_dREDUCEDFORM;
    % check rank criteria for full Jacobian
    if rank_dREDUCEDFORM == totparam_nbr
        % all parameters are identifiable
        no_identification_reducedform = 1; %skip in the following
        indparam_dREDUCEDFORM = [];
    else
        % there is lack of identification
        indparam_dREDUCEDFORM = indtotparam; %initialize for nchoosek
    end
else
    indparam_dREDUCEDFORM = []; %empty for nchoosek
end

% initialize for moments criteria
if ~no_identification_moments && ~error_indicator.identification_moments
    dMOMENTS = ide_moments.dMOMENTS;
    dMOMENTS(ide_moments.ind_dMOMENTS,:) = dMOMENTS(ide_moments.ind_dMOMENTS,:)./ide_moments.norm_dMOMENTS; %normalize
    if strcmp(tol_rank,'robust')
        rank_dMOMENTS = rank(full(dMOMENTS)); %compute rank with imposed tolerance level
    else
        rank_dMOMENTS = rank(full(dMOMENTS),tol_rank); %compute rank with imposed tolerance level
    end
    ide_moments.rank = rank_dMOMENTS;
    % check rank criteria for full Jacobian
    if rank_dMOMENTS == totparam_nbr
        % all parameters are identifiable
        no_identification_moments = 1; %skip in the following
        indparam_dMOMENTS = [];
    else
        % there is lack of identification
        indparam_dMOMENTS = indtotparam; %initialize for nchoosek
    end
else
    indparam_dMOMENTS = []; %empty for nchoosek
end

% initialize for spectrum criteria
if ~no_identification_spectrum && ~error_indicator.identification_spectrum
    dSPECTRUM = ide_spectrum.tilda_dSPECTRUM; %tilda dSPECTRUM is normalized dSPECTRUM matrix in identification.analysis.m
    %alternative normalization
    %dSPECTRUM = ide_spectrum.dSPECTRUM;
    %dSPECTRUM(ide_spectrum.ind_dSPECTRUM,:) = dSPECTRUM(ide_spectrum.ind_dSPECTRUM,:)./ide_spectrum.norm_dSPECTRUM; %normalize
    if strcmp(tol_rank,'robust')
        rank_dSPECTRUM = rank(full(dSPECTRUM)); %compute rank with imposed tolerance level
    else
        rank_dSPECTRUM = rank(full(dSPECTRUM),tol_rank); %compute rank with imposed tolerance level
    end
    ide_spectrum.rank = rank_dSPECTRUM;
    % check rank criteria for full Jacobian
    if rank_dSPECTRUM == totparam_nbr
        % all parameters are identifiable
        no_identification_spectrum = 1; %skip in the following
        indparam_dSPECTRUM = [];
    else
        % lack of identification
        indparam_dSPECTRUM = indtotparam; %initialize for nchoosek
    end
else
    indparam_dSPECTRUM = []; %empty for nchoosek
end

% initialize for minimal system criteria
if ~no_identification_minimal && ~error_indicator.identification_minimal
    dMINIMAL = ide_minimal.dMINIMAL;
    dMINIMAL(ide_minimal.ind_dMINIMAL,:) = dMINIMAL(ide_minimal.ind_dMINIMAL,:)./ide_minimal.norm_dMINIMAL; %normalize
    dMINIMAL_par  = dMINIMAL(:,1:totparam_nbr);       %part of dMINIMAL that is dependent on parameters
    dMINIMAL_rest = dMINIMAL(:,(totparam_nbr+1):end); %part of dMINIMAL that is independent of parameters
    if strcmp(tol_rank,'robust')
        rank_dMINIMAL = rank(full(dMINIMAL)); %compute rank via SVD see function below
    else
        rank_dMINIMAL = rank(full(dMINIMAL),tol_rank); %compute rank via SVD see function below
    end
    ide_minimal.rank = rank_dMINIMAL;
    dMINIMAL_fixed_rank_nbr = size(dMINIMAL_rest,2);
    % check rank criteria for full Jacobian
    if rank_dMINIMAL == totparam_nbr + dMINIMAL_fixed_rank_nbr
        % all parameters are identifiable
        no_identification_minimal = 1; %skip in the following
        indparam_dMINIMAL = [];
    else
        % lack of identification
        indparam_dMINIMAL = indtotparam; %initialize for nchoosek
    end
else
    indparam_dMINIMAL = []; %empty for nchoosek
end

%% Check single parameters
for j=1:totparam_nbr
    if ~no_identification_dynamic
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if strcmp(tol_rank,'robust')
            if rank(full(dDYNAMIC(:,j))) == 0
                dynamic_problpars{1} = [dynamic_problpars{1};j];
            end
        else
            if rank(full(dDYNAMIC(:,j)),tol_rank) == 0
                dynamic_problpars{1} = [dynamic_problpars{1};j];
            end
        end
    end
    if ~no_identification_reducedform && ~error_indicator.identification_reducedform
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if strcmp(tol_rank,'robust')
            if rank(full(dREDUCEDFORM(:,j))) == 0
                reducedform_problpars{1} = [reducedform_problpars{1};j];
            end
        else
            if rank(full(dREDUCEDFORM(:,j)),tol_rank) == 0
                reducedform_problpars{1} = [reducedform_problpars{1};j];
            end
        end
    end
    if ~no_identification_moments && ~error_indicator.identification_moments
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if strcmp(tol_rank,'robust')
            if rank(full(dMOMENTS(:,j))) == 0
                moments_problpars{1} = [moments_problpars{1};j];
            end
        else
            if rank(full(dMOMENTS(:,j)),tol_rank) == 0
                moments_problpars{1} = [moments_problpars{1};j];
            end
        end
    end
    if ~no_identification_spectrum && ~error_indicator.identification_spectrum
        % Diagonal values correspond to single parameters, absolute value needs to be greater than tolerance level
        if abs(dSPECTRUM(j,j)) < tol_rank
            spectrum_problpars{1} = [spectrum_problpars{1};j];
        end
    end
    if ~no_identification_minimal && ~error_indicator.identification_minimal
        % Columns of dMINIMAL_par correspond to single parameters, needs to be augmented with dMINIMAL_rest (part that is independent of parameters),
        % full rank would be equal to 1+dMINIMAL_fixed_rank_nbr
        if strcmp(tol_rank,'robust')
            if rank(full([dMINIMAL_par(:,j) dMINIMAL_rest])) == dMINIMAL_fixed_rank_nbr
                minimal_problpars{1} = [minimal_problpars{1};j];
            end
        else
            if rank(full([dMINIMAL_par(:,j) dMINIMAL_rest]),tol_rank) == dMINIMAL_fixed_rank_nbr
                minimal_problpars{1} = [minimal_problpars{1};j];
            end
        end
    end
end

% Check whether lack of identification is only due to single parameters
if ~no_identification_dynamic
    if size(dynamic_problpars{1},1) == (totparam_nbr - rank_dDYNAMIC)
        %found all nonidentified parameter sets
        no_identification_dynamic = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dDYNAMIC(dynamic_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_reducedform && ~error_indicator.identification_reducedform
    if size(reducedform_problpars{1},1) == (totparam_nbr - rank_dREDUCEDFORM)
        %found all nonidentified parameter sets
        no_identification_reducedform = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dREDUCEDFORM(reducedform_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_moments && ~error_indicator.identification_moments
    if size(moments_problpars{1},1) == (totparam_nbr - rank_dMOMENTS)
        %found all nonidentified parameter sets
        no_identification_moments = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dMOMENTS(moments_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_spectrum && ~error_indicator.identification_spectrum
    if size(spectrum_problpars{1},1) == (totparam_nbr - rank_dSPECTRUM)
        %found all nonidentified parameter sets
        no_identification_spectrum = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dSPECTRUM(spectrum_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_minimal && ~error_indicator.identification_minimal
    if size(minimal_problpars{1},1) == (totparam_nbr + dMINIMAL_fixed_rank_nbr - rank_dMINIMAL)
        %found all nonidentified parameter sets
        no_identification_minimal = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dMINIMAL(minimal_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end


%% check higher order (j>1) parameter sets
%get common parameter indices from which to sample higher-order sets using nchoosek (we do not want to run nchoosek three times), most of the times indparamdMOMENTS, indparamdSPECTRUM, and indparamdMINIMAL are equal anyways
indtotparam = unique([indparam_dDYNAMIC indparam_dREDUCEDFORM indparam_dMOMENTS indparam_dSPECTRUM indparam_dMINIMAL]);

for j=2:min(length(indtotparam),max_dim_subsets_groups) % Check j-element subsets
    h = dyn_waitbar(0,['Brute force collinearity for ' int2str(j) ' parameters.']);
    %Step1: get all possible unique subsets of j elements
    if ~no_identification_dynamic ...
            || (~no_identification_reducedform && ~error_indicator.identification_reducedform)...
            || (~no_identification_moments && ~error_indicator.identification_moments)...
            || (~no_identification_spectrum && ~error_indicator.identification_spectrum)...
            || (~no_identification_minimal && ~error_indicator.identification_minimal)
        indexj=nchoosek(int16(indtotparam),j);  %  int16 speeds up nchoosek
        % One could also use a mex version of nchoosek to speed this up, e.g.VChooseK from https://de.mathworks.com/matlabcentral/fileexchange/26190-vchoosek
    end
    
    %Step 2: remove already problematic sets and initialize rank vector
    if ~no_identification_dynamic
        indexj_dDYNAMIC = RemoveProblematicParameterSets(indexj,dynamic_problpars);
        rankj_dDYNAMIC  = zeros(size(indexj_dDYNAMIC,1),1);
    else
        indexj_dDYNAMIC = [];
    end
    if ~no_identification_reducedform && ~error_indicator.identification_reducedform
        indexj_dREDUCEDFORM = RemoveProblematicParameterSets(indexj,reducedform_problpars);
        rankj_dREDUCEDFORM  = zeros(size(indexj_dREDUCEDFORM,1),1);
    else
        indexj_dREDUCEDFORM = [];
    end
    if ~no_identification_moments && ~error_indicator.identification_moments
        indexj_dMOMENTS = RemoveProblematicParameterSets(indexj,moments_problpars);
        rankj_dMOMENTS = zeros(size(indexj_dMOMENTS,1),1);
    else
        indexj_dMOMENTS = [];
    end
    if ~no_identification_spectrum && ~error_indicator.identification_spectrum
        indexj_dSPECTRUM = RemoveProblematicParameterSets(indexj,spectrum_problpars);
        rankj_dSPECTRUM = zeros(size(indexj_dSPECTRUM,1),1);
    else
        indexj_dSPECTRUM = [];
    end
    if ~no_identification_minimal && ~error_indicator.identification_minimal
        indexj_dMINIMAL = RemoveProblematicParameterSets(indexj,minimal_problpars);
        rankj_dMINIMAL = zeros(size(indexj_dMINIMAL,1),1);
    else
        indexj_dMINIMAL = [];
    end

    %Step3: Check rank criteria on submatrices
    k_dDYNAMIC=0; k_dREDUCEDFORM=0; k_dMOMENTS=0; k_dSPECTRUM=0; k_dMINIMAL=0; %initialize counters
    maxk = max([size(indexj_dDYNAMIC,1), size(indexj_dREDUCEDFORM,1), size(indexj_dMOMENTS,1), size(indexj_dMINIMAL,1), size(indexj_dSPECTRUM,1)]);
    for k=1:maxk
        if ~no_identification_dynamic
            k_dDYNAMIC = k_dDYNAMIC+1;
            if k_dDYNAMIC <= size(indexj_dDYNAMIC,1)
                dDYNAMIC_j = dDYNAMIC(:,indexj_dDYNAMIC(k_dDYNAMIC,:)); % pick columns that correspond to parameter subset
                if strcmp(tol_rank,'robust')
                    rankj_dDYNAMIC(k_dDYNAMIC,1) = rank(full(dDYNAMIC_j)); %compute rank with imposed tolerance
                else
                    rankj_dDYNAMIC(k_dDYNAMIC,1) = rank(full(dDYNAMIC_j),tol_rank); %compute rank with imposed tolerance
                end
            end
        end
        if ~no_identification_reducedform && ~error_indicator.identification_reducedform
            k_dREDUCEDFORM = k_dREDUCEDFORM+1;
            if k_dREDUCEDFORM <= size(indexj_dREDUCEDFORM,1)
                dREDUCEDFORM_j = dREDUCEDFORM(:,indexj_dREDUCEDFORM(k_dREDUCEDFORM,:)); % pick columns that correspond to parameter subset
                if strcmp(tol_rank,'robust')
                    rankj_dREDUCEDFORM(k_dREDUCEDFORM,1) = rank(full(dREDUCEDFORM_j)); %compute rank with imposed tolerance
                else
                    rankj_dREDUCEDFORM(k_dREDUCEDFORM,1) = rank(full(dREDUCEDFORM_j),tol_rank); %compute rank with imposed tolerance
                end
            end
        end
        if ~no_identification_moments && ~error_indicator.identification_moments
            k_dMOMENTS = k_dMOMENTS+1;
            if k_dMOMENTS <= size(indexj_dMOMENTS,1)
                dMOMENTS_j = dMOMENTS(:,indexj_dMOMENTS(k_dMOMENTS,:)); % pick columns that correspond to parameter subset
                if strcmp(tol_rank,'robust')
                    rankj_dMOMENTS(k_dMOMENTS,1) = rank(full(dMOMENTS_j)); %compute rank with imposed tolerance
                else
                    rankj_dMOMENTS(k_dMOMENTS,1) = rank(full(dMOMENTS_j),tol_rank); %compute rank with imposed tolerance
                end
            end
        end
        if ~no_identification_minimal && ~error_indicator.identification_minimal
            k_dMINIMAL = k_dMINIMAL+1;
            if k_dMINIMAL <= size(indexj_dMINIMAL,1)
                dMINIMAL_j = [dMINIMAL_par(:,indexj_dMINIMAL(k_dMINIMAL,:)) dMINIMAL_rest]; % pick columns in dMINIMAL_par that correspond to parameter subset and augment with parameter-indepdendet part dMINIMAL_rest
                if strcmp(tol_rank,'robust')
                    rankj_dMINIMAL(k_dMINIMAL,1) = rank(full(dMINIMAL_j)); %compute rank with imposed tolerance
                else
                    rankj_dMINIMAL(k_dMINIMAL,1) = rank(full(dMINIMAL_j),tol_rank); %compute rank with imposed tolerance
                end
            end
        end
        if ~no_identification_spectrum && ~error_indicator.identification_spectrum
            k_dSPECTRUM = k_dSPECTRUM+1;
            if k_dSPECTRUM <= size(indexj_dSPECTRUM,1)
                dSPECTRUM_j = dSPECTRUM(indexj_dSPECTRUM(k_dSPECTRUM,:),indexj_dSPECTRUM(k_dSPECTRUM,:)); % pick rows and columns that correspond to parameter subset
                if strcmp(tol_rank,'robust')
                    rankj_dSPECTRUM(k_dSPECTRUM,1) = rank(full(dSPECTRUM_j)); % Compute rank with imposed tol
                else
                    rankj_dSPECTRUM(k_dSPECTRUM,1) = rank(full(dSPECTRUM_j),tol_rank); % Compute rank with imposed tol
                end
            end
        end
        dyn_waitbar(k/maxk,h)
    end

    %Step 4: Compare rank conditions for all possible subsets. If rank condition is violated, then the corresponding numbers of the parameters are stored
    if ~no_identification_dynamic
        dynamic_problpars{j} = indexj_dDYNAMIC(rankj_dDYNAMIC < j,:);
    end
    if ~no_identification_reducedform && ~error_indicator.identification_reducedform
        reducedform_problpars{j} = indexj_dREDUCEDFORM(rankj_dREDUCEDFORM < j,:);
    end
    if ~no_identification_moments && ~error_indicator.identification_moments
        moments_problpars{j} = indexj_dMOMENTS(rankj_dMOMENTS < j,:);
    end
    if ~no_identification_minimal && ~error_indicator.identification_minimal
        minimal_problpars{j} = indexj_dMINIMAL(rankj_dMINIMAL < (j+dMINIMAL_fixed_rank_nbr),:);
    end
    if ~no_identification_spectrum && ~error_indicator.identification_spectrum
        spectrum_problpars{j} = indexj_dSPECTRUM(rankj_dSPECTRUM < j,:);
    end
%     % Optional Step 5: % remove redundant 2-sets, eg. if the problematic sets are [(p1,p2);(p1,p3);(p2,p3)], then the unique problematic parameter sets are actually only [(p1,p2),(p1,p3)]        
%     if j == 2        
%         for jj=1:max([size(dynamic_problpars{2},1), size(reducedform_problpars{2},1), size(moments_problpars{2},1), size(spectrum_problpars{2},1), size(minimal_problpars{2},1)])
%             if jj <= size(dynamic_problpars{2},1)
%                 dynamic_problpars{2}(dynamic_problpars{2}(jj,2)==dynamic_problpars{2}(:,1)) = dynamic_problpars{2}(jj,1);
%             end
%             if jj <= size(reducedform_problpars{2},1)
%                 reducedform_problpars{2}(reducedform_problpars{2}(jj,2)==reducedform_problpars{2}(:,1)) = reducedform_problpars{2}(jj,1);
%             end
%             if jj <= size(moments_problpars{2},1)
%                 moments_problpars{2}(moments_problpars{2}(jj,2)==moments_problpars{2}(:,1)) = moments_problpars{2}(jj,1);
%             end
%             if jj <= size(spectrum_problpars{2},1)
%                 spectrum_problpars{2}(spectrum_problpars{2}(jj,2)==spectrum_problpars{2}(:,1)) = spectrum_problpars{2}(jj,1);
%             end
%             if jj <= size(minimal_problpars{2},1)
%                 minimal_problpars{2}(minimal_problpars{2}(jj,2)==minimal_problpars{2}(:,1)) = minimal_problpars{2}(jj,1);
%             end
%         end
%         dynamic_problpars{2} = unique(dynamic_problpars{2},'rows');
%         reducedform_problpars{2} = unique(reducedform_problpars{2},'rows');
%         moments_problpars{2} = unique(moments_problpars{2},'rows');
%         spectrum_problpars{2} = unique(spectrum_problpars{2},'rows');
%         minimal_problpars{2} = unique(minimal_problpars{2},'rows');
%         % in indparam we replace the second parameter of problematic 2-sets by the observational equivalent first parameter to speed up nchoosek
%         idx2 = unique([dynamic_problpars{2}; reducedform_problpars{2}; moments_problpars{2}; spectrum_problpars{2}; minimal_problpars{2}],'rows');
%         if ~isempty(idx2)
%             indtotparam(ismember(indtotparam,idx2(:,2))) = [];
%         end
%     end
    dyn_waitbar_close(h);
end

%% Save output variables
if ~isempty(dynamic_problpars{1})
    dynamic_problpars{1}(ismember(dynamic_problpars{1},1:(totparam_nbr-modparam_nbr))) = []; % get rid of stderr and corr variables for dynamic
end
ide_dynamic.problpars         = dynamic_problpars;
ide_reducedform.problpars = reducedform_problpars;
ide_moments.problpars     = moments_problpars;
ide_spectrum.problpars    = spectrum_problpars;
ide_minimal.problpars     = minimal_problpars;

%% Auxiliary functions
function idx = RemoveProblematicParameterSets(idx,problparset)
% Remove already problematic parameters
% INPUTS:
%   * idx:         complete index of possible combinations
%   * problparset  [cell] of all lower order combinations that are problematic
%   * iset         [integer] number of elements in set to consider
% OUTPUTS:
%   * idx:         index of possible combinations without already
%                  problematic lower order sets
    iset = size(idx,2);
    for iii=1:(iset-1)
        if ~isempty(problparset{iii})
            for kkk=1:size(problparset{iii},1)                
                idx((sum(ismember(idx,problparset{iii}(kkk,:)),2)==iii),:) = [];
            end
        end
    end 
end%RemoveProblematicParameterSets ed


end %main function end