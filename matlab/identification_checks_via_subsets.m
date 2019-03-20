function [ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = identification_checks_via_subsets(ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident)
%[ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = identification_checks_via_subsets(ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident)
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
%                       and Iskrev, 2011). If ide_reducedform.no_identification_reducedform
%                       is 1 then the search for problematic parameter sets will be skipped
%   ide_moments:        [structure] Containing results from identification
%                       analysis based on moments (Iskrev, 2010). If
%                       ide_moments.no_identification_moments is 1 then the search for
%                       problematic parameter sets will be skipped
%   ide_spectrum:       [structure] Containing results from identification
%                       analysis based on the spectrum (Qu and Tkachenko, 2012).
%                       If ide_spectrum.no_identification_spectrum is 1 then the search for
%                       problematic parameter sets will be skipped
%   ide_minimal:        [structure] Containing results from identification
%                       analysis based on the minimal state space system
%                       (Komunjer and Ng, 2011). If ide_minimal.no_identification_minimal
%                       is 1 then the search for problematic parameter sets will be skipped
%   totparam_nbr:       [integer] number of estimated stderr, corr and model parameters
%   numzerotolrank:     [double] tolerance level for rank compuations
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
%   * identification_analysis.m 
% =========================================================================
% Copyright (C) 2019 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

%% initialize output objects and get options
no_identification_lre                 = 0; %always compute lre
no_identification_reducedform         = options_ident.no_identification_reducedform;
no_identification_moments             = options_ident.no_identification_moments;
no_identification_spectrum            = options_ident.no_identification_spectrum;
no_identification_minimal             = options_ident.no_identification_minimal;
tol_rank               = options_ident.tol_rank;
max_dim_subsets_groups = options_ident.max_dim_subsets_groups;
lre_problpars          = cell(1,max_dim_subsets_groups);
reducedform_problpars  = cell(1,max_dim_subsets_groups);
moments_problpars      = cell(1,max_dim_subsets_groups);
spectrum_problpars     = cell(1,max_dim_subsets_groups);
minimal_problpars      = cell(1,max_dim_subsets_groups);

indtotparam            = 1:totparam_nbr; %initialize index of parameters

%% Prepare Jacobians and check rank
% initialize linear rational expectations model
if ~no_identification_lre
    dLRE = ide_lre.dLRE;
    dLRE(ide_lre.ind_dLRE,:) = dLRE(ide_lre.ind_dLRE,:)./ide_lre.norm_dLRE; %normalize
    if totparam_nbr > modparam_nbr
        dLRE = [zeros(size(ide_lre.dLRE,1),totparam_nbr-modparam_nbr) dLRE]; %add derivatives wrt stderr and corr parameters
    end
    rank_dLRE = rank(dLRE,tol_rank); %compute rank with imposed tolerance level
    ide_lre.rank = rank_dLRE;
    % check rank criteria for full Jacobian
    if rank_dLRE == totparam_nbr 
        % all parameters are identifiable
        no_identification_lre = 1; %skip in the following
        indparam_dLRE = [];
    else 
        % there is lack of identification
        indparam_dLRE = indtotparam; %initialize for nchoosek
    end
else
    indparam_dLRE = []; %empty for nchoosek
end

% initialize for reduced form solution criteria
if ~no_identification_reducedform
    dTAU = ide_reducedform.dTAU;
    dTAU(ide_reducedform.ind_dTAU,:) = dTAU(ide_reducedform.ind_dTAU,:)./ide_reducedform.norm_dTAU; %normalize
    rank_dTAU = rank(dTAU,tol_rank); %compute rank with imposed tolerance level
    ide_reducedform.rank = rank_dTAU;
    % check rank criteria for full Jacobian
    if rank_dTAU == totparam_nbr
        % all parameters are identifiable
        no_identification_reducedform = 1; %skip in the following
        indparam_dTAU = [];
    else 
        % there is lack of identification
        indparam_dTAU = indtotparam; %initialize for nchoosek
    end
else
    indparam_dTAU = []; %empty for nchoosek
end

% initialize for moments criteria
if ~no_identification_moments
    J = ide_moments.J;
    J(ide_moments.ind_J,:) = J(ide_moments.ind_J,:)./ide_moments.norm_J; %normalize
    rank_J = rank(J,tol_rank); %compute rank with imposed tolerance level
    ide_moments.rank = rank_J;
    % check rank criteria for full Jacobian
    if rank_J == totparam_nbr 
        % all parameters are identifiable
        no_identification_moments = 1; %skip in the following
        indparam_J = [];
    else 
        % there is lack of identification
        indparam_J = indtotparam; %initialize for nchoosek
    end
else
    indparam_J = []; %empty for nchoosek
end

% initialize for spectrum criteria
if ~no_identification_spectrum
    G = ide_spectrum.tilda_G; %tilda G is normalized G matrix in identification_analysis.m
    %alternative normalization
    %G = ide_spectrum.G;
    %G(ide_spectrum.ind_G,:) = G(ide_spectrum.ind_G,:)./ide_spectrum.norm_G; %normalize
    rank_G = rank(G,tol_rank); %compute rank with imposed tolerance level
    ide_spectrum.rank = rank_G;
    % check rank criteria for full Jacobian
    if rank_G == totparam_nbr
        % all parameters are identifiable
        no_identification_spectrum = 1; %skip in the following
        indparam_G = [];
    else
        % lack of identification
        indparam_G = indtotparam; %initialize for nchoosek
    end
else
    indparam_G = []; %empty for nchoosek
end

% initialize for minimal system criteria
if ~no_identification_minimal
    D = ide_minimal.D;
    D(ide_minimal.ind_D,:) = D(ide_minimal.ind_D,:)./ide_minimal.norm_D; %normalize
    D_par  = D(:,1:totparam_nbr);       %part of D that is dependent on parameters
    D_rest = D(:,(totparam_nbr+1):end); %part of D that is independent of parameters 
    rank_D = rank(D,tol_rank); %compute rank via SVD see function below
    ide_minimal.rank = rank_D;
    D_fixed_rank_nbr = size(D_rest,2);
    % check rank criteria for full Jacobian
    if rank_D == totparam_nbr + D_fixed_rank_nbr
        % all parameters are identifiable
        no_identification_minimal = 1; %skip in the following
        indparam_D = [];
    else
        % lack of identification
        indparam_D = indtotparam; %initialize for nchoosek
    end
else
    indparam_D = []; %empty for nchoosek
end

%% Check single parameters
for j=1:totparam_nbr
    if ~no_identification_lre
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if rank(dLRE(:,j),tol_rank) == 0
            lre_problpars{1} = [lre_problpars{1};j];
        end
    end
    if ~no_identification_reducedform
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if rank(dTAU(:,j),tol_rank) == 0
            reducedform_problpars{1} = [reducedform_problpars{1};j];
        end
    end
    if ~no_identification_moments
        % Columns correspond to single parameters, i.e. full rank would be equal to 1
        if rank(J(:,j),tol_rank) == 0
            moments_problpars{1} = [moments_problpars{1};j];
        end
    end
    if ~no_identification_spectrum
        % Diagonal values correspond to single parameters, absolute value needs to be greater than tolerance level
        if abs(G(j,j)) < tol_rank
            spectrum_problpars{1} = [spectrum_problpars{1};j];
        end
    end
    if ~no_identification_minimal
        % Columns of D_par correspond to single parameters, needs to be augmented with D_rest (part that is independent of parameters),
        % full rank would be equal to 1+D_fixed_rank_nbr
        if rank([D_par(:,j) D_rest],tol_rank) == D_fixed_rank_nbr
            minimal_problpars{1} = [minimal_problpars{1};j];
        end
    end
end

% Check whether lack of identification is only due to single parameters
if ~no_identification_lre
    if size(lre_problpars{1},1) == (totparam_nbr - rank_dLRE)
        %found all nonidentified parameter sets
        no_identification_lre = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dLRE(lre_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_reducedform
    if size(reducedform_problpars{1},1) == (totparam_nbr - rank_dTAU)
        %found all nonidentified parameter sets
        no_identification_reducedform = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_dTAU(reducedform_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_moments
    if size(moments_problpars{1},1) == (totparam_nbr - rank_J)
        %found all nonidentified parameter sets
        no_identification_moments = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_J(moments_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_spectrum
    if size(spectrum_problpars{1},1) == (totparam_nbr - rank_G)
        %found all nonidentified parameter sets
        no_identification_spectrum = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_G(spectrum_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end
if ~no_identification_minimal
    if size(minimal_problpars{1},1) == (totparam_nbr + D_fixed_rank_nbr - rank_D)
        %found all nonidentified parameter sets
        no_identification_minimal = 1; %skip in the following
    else
        %still parameter sets that are nonidentified
        indparam_D(minimal_problpars{1}) = []; %remove single unidentified parameters from higher-order sets of indparam
    end
end


%% check higher order (j>1) parameter sets
%get common parameter indices from which to sample higher-order sets using nchoosek (we do not want to run nchoosek three times), most of the times indparamJ, indparamG, and indparamD are equal anyways
indtotparam = unique([indparam_dLRE indparam_dTAU indparam_J indparam_G indparam_D]);

for j=2:min(length(indtotparam),max_dim_subsets_groups) % Check j-element subsets
    h = dyn_waitbar(0,['Brute force collinearity for ' int2str(j) ' parameters.']);
    %Step1: get all possible unique subsets of j elements
    if ~no_identification_lre || ~no_identification_reducedform || ~no_identification_moments || ~no_identification_spectrum || ~no_identification_minimal
        indexj=nchoosek(int16(indtotparam),j);  %  int16 speeds up nchoosek
        % One could also use a mex version of nchoosek to speed this up, e.g.VChooseK from https://de.mathworks.com/matlabcentral/fileexchange/26190-vchoosek
    end
    
    %Step 2: remove already problematic sets and initialize rank vector
    if ~no_identification_lre
        indexj_dLRE = RemoveProblematicParameterSets(indexj,lre_problpars);
        rankj_dLRE  = zeros(size(indexj_dLRE,1),1);
    else
        indexj_dLRE = [];
    end
    if ~no_identification_reducedform
        indexj_dTAU = RemoveProblematicParameterSets(indexj,reducedform_problpars);
        rankj_dTAU  = zeros(size(indexj_dTAU,1),1);
    else
        indexj_dTAU = [];
    end
    if ~no_identification_moments
        indexj_J = RemoveProblematicParameterSets(indexj,moments_problpars);
        rankj_J = zeros(size(indexj_J,1),1);
    else
        indexj_J = [];
    end
    if ~no_identification_spectrum
        indexj_G = RemoveProblematicParameterSets(indexj,spectrum_problpars);
        rankj_G = zeros(size(indexj_G,1),1);
    else
        indexj_G = [];
    end
    if ~no_identification_minimal
        indexj_D = RemoveProblematicParameterSets(indexj,minimal_problpars);
        rankj_D = zeros(size(indexj_D,1),1);
    else
        indexj_D = [];
    end
    
    %Step3: Check rank criteria on submatrices
    k_dLRE=0; k_dTAU=0; k_J=0; k_G=0; k_D=0; %initialize counters
    maxk = max([size(indexj_dLRE,1), size(indexj_dTAU,1), size(indexj_J,1), size(indexj_D,1), size(indexj_G,1)]);
    for k=1:maxk
        if ~no_identification_lre
            k_dLRE = k_dLRE+1;
            if k_dLRE <= size(indexj_dLRE,1)
                dLRE_j = dLRE(:,indexj_dLRE(k_dLRE,:)); % pick columns that correspond to parameter subset
                rankj_dLRE(k_dLRE,1) = rank(dLRE_j,tol_rank); %compute rank with imposed tolerance
            end
        end
        if ~no_identification_reducedform
            k_dTAU = k_dTAU+1;
            if k_dTAU <= size(indexj_dTAU,1)
                dTAU_j = dTAU(:,indexj_dTAU(k_dTAU,:)); % pick columns that correspond to parameter subset
                rankj_dTAU(k_dTAU,1) = rank(dTAU_j,tol_rank); %compute rank with imposed tolerance
            end
        end
        if ~no_identification_moments
            k_J = k_J+1;
            if k_J <= size(indexj_J,1)
                J_j = J(:,indexj_J(k_J,:)); % pick columns that correspond to parameter subset
                rankj_J(k_J,1) = rank(J_j,tol_rank); %compute rank with imposed tolerance
            end
        end
        if ~no_identification_minimal
            k_D = k_D+1;
            if k_D <= size(indexj_D,1)
                D_j = [D_par(:,indexj_D(k_D,:)) D_rest]; % pick columns in D_par that correspond to parameter subset and augment with parameter-indepdendet part D_rest
                rankj_D(k_D,1) = rank(D_j,tol_rank); %compute rank with imposed tolerance
            end
        end
        if ~no_identification_spectrum
            k_G = k_G+1;
            if k_G <= size(indexj_G,1)
                G_j = G(indexj_G(k_G,:),indexj_G(k_G,:)); % pick rows and columns that correspond to parameter subset
                rankj_G(k_G,1) = rank(G_j,tol_rank); % Compute rank with imposed tol
            end
        end
        dyn_waitbar(k/maxk,h)
    end
    
    %Step 4: Compare rank conditions for all possible subsets. If rank condition is violated, then the corresponding numbers of the parameters are stored
    if ~no_identification_lre
        lre_problpars{j} = indexj_dLRE(rankj_dLRE < j,:);
    end
    if ~no_identification_reducedform
        reducedform_problpars{j} = indexj_dTAU(rankj_dTAU < j,:);
    end
    if ~no_identification_moments        
        moments_problpars{j} = indexj_J(rankj_J < j,:);
    end
    if ~no_identification_minimal
        minimal_problpars{j} = indexj_D(rankj_D < (j+D_fixed_rank_nbr),:);
    end
    if ~no_identification_spectrum
        spectrum_problpars{j} = indexj_G(rankj_G < j,:);
    end
%     % Optional Step 5: % remove redundant 2-sets, eg. if the problematic sets are [(p1,p2);(p1,p3);(p2,p3)], then the unique problematic parameter sets are actually only [(p1,p2),(p1,p3)]        
%     if j == 2        
%         for jj=1:max([size(lre_problpars{2},1), size(reducedform_problpars{2},1), size(moments_problpars{2},1), size(spectrum_problpars{2},1), size(minimal_problpars{2},1)])
%             if jj <= size(lre_problpars{2},1)
%                 lre_problpars{2}(lre_problpars{2}(jj,2)==lre_problpars{2}(:,1)) = lre_problpars{2}(jj,1);
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
%         lre_problpars{2} = unique(lre_problpars{2},'rows');
%         reducedform_problpars{2} = unique(reducedform_problpars{2},'rows');
%         moments_problpars{2} = unique(moments_problpars{2},'rows');
%         spectrum_problpars{2} = unique(spectrum_problpars{2},'rows');
%         minimal_problpars{2} = unique(minimal_problpars{2},'rows');
%         % in indparam we replace the second parameter of problematic 2-sets by the observational equivalent first parameter to speed up nchoosek
%         idx2 = unique([lre_problpars{2}; reducedform_problpars{2}; moments_problpars{2}; spectrum_problpars{2}; minimal_problpars{2}],'rows');
%         if ~isempty(idx2)
%             indtotparam(ismember(indtotparam,idx2(:,2))) = [];
%         end
%     end
    dyn_waitbar_close(h);
end

%% Save output variables
if ~isempty(lre_problpars{1})
    lre_problpars{1}(ismember(lre_problpars{1},1:(totparam_nbr-modparam_nbr))) = []; % get rid of stderr and corr variables for lre
end
ide_lre.problpars         = lre_problpars;
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