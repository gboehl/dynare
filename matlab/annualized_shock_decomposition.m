function [z, endo_names, endo_names_tex, steady_state, i_var, oo_] = annualized_shock_decomposition(oo_, M_, opts, i_var, t0, t1, realtime_, vintage_, steady_state,GYTREND0,var_type,islog)
% function oo_ = annualized_shock_decomposition(oo_,t0,options_.nobs);
% Computes annualized shocks contribution to a simulated trajectory. The fields set are
% oo_.annualized_shock_decomposition, oo_.annualized_realtime_shock_decomposition, 
% oo_.annualized_realtime_conditional_shock_decomposition and oo_.annualized_realtime_forecast_shock_decomposition. 
% Subfields are arrays n_var by nshock+2 by nperiods. The
% first nshock columns store the respective shock contributions, column n+1
% stores the role of the initial conditions, while column n+2 stores the
% value of the smoothed variables.  Both the variables and shocks are stored 
% in the order of endo_names and M_.exo_names, respectively.
%
% INPUTS
%    oo_:          [structure] Storage of results
%    M_:           [structure] Storage of model
%    opts:         [structure] options for shock decomp 
%    i_var:        [array] index of vars
%    t0:           [integer]  first period
%    t1:           [integer] last period
%    realtime_:    [integer]   
%    vintage_:     [integer]
%    steady_state: [array] steady state value of quarterly (log-) level vars
%    GYTREND0:     [array] growth of level of vars
%    var_type:     [integer] flag for stock/flow/deflator
%    islog:        [integer] flag for log-levels 
%
% OUTPUTS
%    z:              [matrix] shock decomp to plot 
%    endo_names:     [char] updated var names
%    endo_names_tex: [char] updated TeX var names
%    steady_state:   [array] updated stady state of vars
%    i_var:          [integer array] updated var indices to plot
%    oo_:            [structure]  Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

nvar = length(i_var);

% usual shock decomp 
z = oo_.shock_decomposition;
z = z(i_var,:,:);
if var_type==2,
    gtxt = 'PHI'; % inflation rate
    gtex = '\pi';
else
    gtxt = 'G'; % growth rate
    gtex = 'g';
end
steady_state=steady_state(i_var);
% endo_names = M_.endo_names(i_var,:);
% endo_names_tex = M_.endo_names_tex(i_var,:);
nterms = size(z,2);
nfrcst = opts.forecast/4;

for j=1:nvar
    if j>1,
        endo_names = char(endo_names,[deblank(M_.endo_names(i_var(j),:)) '_A']);
        endo_names_tex = char(endo_names_tex,['{' deblank(M_.endo_names_tex(i_var(j),:)) '}^A']);
        gendo_names = char(gendo_names,[gtxt endo_names(j,:)]);
        gendo_names_tex = char(gendo_names_tex,[gtex '(' deblank(endo_names_tex(j,:)) ')']);
    else
        endo_names = [deblank(M_.endo_names(i_var(j),:)) '_A'];
        endo_names_tex = ['{' deblank(M_.endo_names_tex(i_var(j),:)) '}^A'];
        gendo_names = [gtxt endo_names(j,:)];
        gendo_names_tex = [gtex '(' deblank(endo_names_tex(j,:)) ')'];
    end
    for k =1:nterms,
        [za(j,k,:), steady_state_a(j,1), gza(j,k,:), steady_state_ga(j,1)] = ...
            quarterly2annual(squeeze(z(j,k,t0:end)),steady_state(j),GYTREND0,var_type,islog);
    end
    ztmp=squeeze(za(j,:,:));
    zscale = sum(ztmp(1:end-1,:))./ztmp(end,:);
    ztmp(1:end-1,:) = ztmp(1:end-1,:)./repmat(zscale,[nterms-1,1]);
    gztmp=squeeze(gza(j,:,:));
    gscale = sum(gztmp(1:end-1,:))./ gztmp(end,:);
    gztmp(1:end-1,:) = gztmp(1:end-1,:)./repmat(gscale,[nterms-1,1]);
    za(j,:,:) = ztmp;
    gza(j,:,:) = gztmp;
end

z=cat(1,za,gza);
oo_.annualized_shock_decomposition=z;
endo_names = char(endo_names,gendo_names);
endo_names_tex = char(endo_names_tex,gendo_names_tex);

% realtime
init=1;
for i=t0:4:t1,
    yr=floor(i/4);
    za=[];
    gza=[];
    z = oo_.realtime_shock_decomposition.(['time_' int2str(i)]);
    z = z(i_var,:,:);
           
    for j=1:nvar
        for k =nterms:-1:1,
%             if k<nterms
%                 ztmp = squeeze(sum(z(j,[1:k-1,k+1:end-1],t0-4:end)));
%             else
                ztmp = squeeze(z(j,k,t0-4:end));
%             end
            [za(j,k,:), steady_state_a(j,1), gza(j,k,:), steady_state_ga(j,1)] = ...
                quarterly2annual(ztmp,steady_state(j),GYTREND0,var_type,islog);
%             if k<nterms
%                 za(j,k,:) = za(j,end,:) - za(j,k,:);
%                 gza(j,k,:) = gza(j,end,:) - gza(j,k,:);
%             end
            
        end
        
        ztmp=squeeze(za(j,:,:));
        gztmp=squeeze(gza(j,:,:));
        ztmp=squeeze(za(j,:,:));
        zscale = sum(ztmp(1:end-1,:))./ztmp(end,:);
        ztmp(1:end-1,:) = ztmp(1:end-1,:)./repmat(zscale,[nterms-1,1]);
        gztmp=squeeze(gza(j,:,:));
    gscale = sum(gztmp(1:end-1,:))./ gztmp(end,:);
        gztmp(1:end-1,:) = gztmp(1:end-1,:)./repmat(gscale,[nterms-1,1]);
        za(j,:,:) = ztmp;
        gza(j,:,:) = gztmp;
    end
    
    z=cat(1,za,gza);
    steady_state = [steady_state_a;steady_state_ga];
    if init==1,
        oo_.annualized_realtime_shock_decomposition.pool = z;
    else
        oo_.annualized_realtime_shock_decomposition.pool(:,:,yr) = z(:,:,end-nfrcst);
    end        
    oo_.annualized_realtime_shock_decomposition.(['yr_' int2str(yr)]) = z;
        
    if opts.forecast
        oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr)]) = z(:,:,end-nfrcst:end);
        if init>nfrcst
            oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)]) = ...
                oo_.annualized_realtime_shock_decomposition.pool(:,:,yr-nfrcst:end) - ...
                oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr-nfrcst)]);
            oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end-1,:) = ...
                oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(yr-nfrcst)])(:,end,:);
        end
    end
% ztmp=oo_.realtime_shock_decomposition.pool(:,:,21:29)-oo_.realtime_forecast_shock_decomposition.time_21;



    init=init+1;
end

i_var=1:2*nvar;
steady_state = [steady_state_a;steady_state_ga];


switch realtime_
    
    case 0
        z = oo_.annualized_shock_decomposition;
    
    case 1 % realtime
        if vintage_
            z = oo_.annualized_realtime_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            z = oo_.annualized_realtime_shock_decomposition.pool;
        end
    
    case 2 % conditional
        if vintage_
            z = oo_.annualized_realtime_conditional_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            error();
        end
        
    case 3 % forecast
        if vintage_
            z = oo_.annualized_realtime_forecast_shock_decomposition.(['yr_' int2str(floor(vintage_/4))]);
        else
            error()
        end
end
