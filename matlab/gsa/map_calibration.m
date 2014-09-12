function map_calibration(OutputDirectoryName, Model, DynareOptions, DynareResults, EstimatedParameters)

% Copyright (C) 2014 Dynare Team
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

fname_ = Model.fname;
pnames = Model.param_names(EstimatedParameters.param_vals(:,1),:);
pvalue_ks = DynareOptions.opt_gsa.pvalue_ks;
indx_irf = [];
indx_moment = [];

skipline()
disp('Sensitivity analysis for calibration criteria')

filetoload=[OutputDirectoryName '/' fname_ '_prior'];
load(filetoload,'lpmat','lpmat0','istable','iunstable','iindeterm','iwrong' ,'infox')
if ~isempty(lpmat0),
    lpmatx=lpmat0(istable,:);
else
    lpmatx=[];
end
[Nsam, np] = size(lpmat);
nshock = size(lpmat0,2);
npT = np+nshock;

nbr_irf_restrictions = size(DynareOptions.endogenous_prior_restrictions.irf,1);
mat_irf=cell(nbr_irf_restrictions,1);
for ij=1:nbr_irf_restrictions,
    mat_irf{ij}=NaN(Nsam,length(DynareOptions.endogenous_prior_restrictions.irf{ij,3}));
end

nbr_moment_restrictions = size(DynareOptions.endogenous_prior_restrictions.moment,1);
mat_moment=cell(nbr_moment_restrictions,1);
for ij=1:nbr_moment_restrictions,
    mat_moment{ij}=NaN(Nsam,length(DynareOptions.endogenous_prior_restrictions.moment{ij,3}));
end

irestrictions = [1:Nsam];
for j=1:Nsam,
    Model = set_all_parameters([lpmat0(j,:) lpmat(j,:)]',EstimatedParameters,Model);
    [Tt,Rr,SteadyState,info] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');
    if info(1)==0,
        [info, info_irf, info_moment, data_irf, data_moment]=endogenous_prior_restrictions(Tt,Rr,Model,DynareOptions,DynareResults);
        if ~isempty(info_irf)
            for ij=1:nbr_irf_restrictions,
                mat_irf{ij}(j,:)=data_irf{ij}(:,2)';
            end
            indx_irf(j,:)=info_irf(:,1);
        end
        if ~isempty(info_moment)
            for ij=1:nbr_moment_restrictions,
                mat_moment{ij}(j,:)=data_moment{ij}(:,2)';
            end
            indx_moment(j,:)=info_moment(:,1);
        end
    else
        irestrictions(j)=0;
    end
end
irestrictions=irestrictions(find(irestrictions));
xmat=[lpmat0(irestrictions,:) lpmat(irestrictions,:)];
skipline()

if ~isempty(indx_irf),
    indx_irf = indx_irf(irestrictions,:);
    h1=dyn_figure(DynareOptions,'name','Evaluation of irf restrictions');
    nrow=ceil(sqrt(nbr_irf_restrictions));
    ncol=nrow;
    if nrow*(nrow-1)>nbr_irf_restrictions,
        ncol=nrow-1;
    end
    for ij=1:nbr_irf_restrictions,
        figure(h1),
        mat_irf{ij}=mat_irf{ij}(irestrictions,:);
        subplot(nrow,ncol,ij),
        for ik=1:size(mat_irf{ij},2),
            cumplot(mat_irf{ij}(:,ik)),
            hold all,
        end
        %     hist(mat_irf{ij}),
        a=axis;
        x1val=max(DynareOptions.endogenous_prior_restrictions.irf{ij,4}(1),a(1));
        x2val=min(DynareOptions.endogenous_prior_restrictions.irf{ij,4}(2),a(2));
        hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'k');
        set(hp,'FaceAlpha',[0.5])
        hold off,
        leg = num2str(DynareOptions.endogenous_prior_restrictions.irf{ij,3}(1));
        if size(mat_irf{ij},2)>1,
            leg = [leg,':' ,num2str(DynareOptions.endogenous_prior_restrictions.irf{ij,3}(end))];
        end
        title([DynareOptions.endogenous_prior_restrictions.irf{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.irf{ij,2}, '(', leg,')'],'interpreter','none'),
        
        indx1 = find(indx_irf(:,ij)==0);
        indx2 = find(indx_irf(:,ij)~=0);
        atitle=[DynareOptions.endogenous_prior_restrictions.irf{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.irf{ij,2}, '(', leg,')'];
        fprintf(['%4.1f%% of the prior support matches IRF ',atitle,' inside [%4.1f, %4.1f]\n'],length(indx1)/length(irestrictions)*100,DynareOptions.endogenous_prior_restrictions.irf{ij,4})
        aname=['irf_calib_',int2str(ij)];
        atitle=['IRF Calib: Parameter(s) driving ',DynareOptions.endogenous_prior_restrictions.irf{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.irf{ij,2}, '(', leg,')'];
        [proba, dproba] = stab_map_1(xmat, indx1, indx2, aname, 0);
        indplot=find(proba<pvalue_ks);
        if ~isempty(indplot)
            stab_map_1(xmat, indx1, indx2, aname, 1, indplot, OutputDirectoryName,[],atitle);
        end
    end
    dyn_saveas(h1,[OutputDirectoryName,filesep,fname_,'_irf_restrictions'],DynareOptions);
    skipline()
end

if ~isempty(indx_moment)
    indx_moment = indx_moment(irestrictions,:);
    h2=dyn_figure(DynareOptions,'name','Evaluation of moment restrictions');
    nrow=ceil(sqrt(nbr_moment_restrictions));
    ncol=nrow;
    if nrow*(nrow-1)>nbr_moment_restrictions,
        ncol=nrow-1;
    end
    for ij=1:nbr_moment_restrictions,
        figure(h2),
        mat_moment{ij}=mat_moment{ij}(irestrictions,:);
        subplot(nrow,ncol,ij),
        for ik=1:size(mat_moment{ij},2),
            cumplot(mat_moment{ij}(:,ik)),
            hold all,
        end
        %     hist(mat_moment{ij}),
        a=axis;
        x1val=max(DynareOptions.endogenous_prior_restrictions.moment{ij,4}(1),a(1));
        x2val=min(DynareOptions.endogenous_prior_restrictions.moment{ij,4}(2),a(2));
        hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'k');
        set(hp,'FaceAlpha',[0.5])
        hold off,
        leg = num2str(DynareOptions.endogenous_prior_restrictions.moment{ij,3}(1));
        if size(mat_moment{ij},2)>1,
            leg = [leg,':' ,num2str(DynareOptions.endogenous_prior_restrictions.moment{ij,3}(end))];
        end
        title([DynareOptions.endogenous_prior_restrictions.moment{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.moment{ij,2},'(',leg,')'],'interpreter','none'),
        
        indx1 = find(indx_moment(:,ij)==0);
        indx2 = find(indx_moment(:,ij)~=0);
        atitle=[DynareOptions.endogenous_prior_restrictions.moment{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.moment{ij,2}, '(', leg,')'];
        fprintf(['%4.1f%% of the prior support matches MOMENT ',atitle,' inside [%4.1f, %4.1f]\n'],length(indx1)/length(irestrictions)*100,DynareOptions.endogenous_prior_restrictions.moment{ij,4})
        aname=['moment_calib_',int2str(ij)];
        atitle=['MOMENT Calib: Parameter(s) driving ',DynareOptions.endogenous_prior_restrictions.moment{ij,1},' vs ',DynareOptions.endogenous_prior_restrictions.moment{ij,2}, '(', leg,')'];
        [proba, dproba] = stab_map_1(xmat, indx1, indx2, aname, 0);
        indplot=find(proba<pvalue_ks);
        if ~isempty(indplot)
            stab_map_1(xmat, indx1, indx2, aname, 1, indplot, OutputDirectoryName,[],atitle);
        end
    end
    dyn_saveas(h2,[OutputDirectoryName,filesep,fname_,'_moment_restrictions'],DynareOptions);
    skipline()
end
return

