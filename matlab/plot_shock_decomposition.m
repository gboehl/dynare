function plot_shock_decomposition(M_,oo_,options_,varlist)
% function plot_shock_decomposition(M_,oo_,options_,varlist)
% Plots the results of shock_decomposition
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2016-2017 Dynare Team
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

% indices of endogenous variables
if size(varlist,1) == 0
    varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
end

[i_var,nvar] = varlist_indices(varlist,M_.endo_names);

% number of variables
endo_nbr = M_.endo_nbr;

% number of shocks
nshocks = M_.exo_nbr;
% type = '';
% fig_names='';
% detail_plot=0;
% realtime_=0; % 0 is standard; 1 is realtime (pool/vintage); 2 is conditional (pool/vintage); 3 is forecast (pool/vintage)
% vintage_=0; % 0 pool realtime/conditional; int: forecast/conditional shock decompositions
% forecast_=0;
% steadystate=0;
% write_xls=0;

    if ~isempty(options_.shock_decomp.fig_names)
        fig_names=['_' options_.shock_decomp.fig_names];
    end
    type=options_.shock_decomp.type;
    detail_plot=options_.shock_decomp.detail_plot;
    realtime_= options_.shock_decomp.realtime;
    vintage_ = options_.shock_decomp.vintage;
    forecast_ = options_.shock_decomp.forecast;
    steadystate = options_.shock_decomp.steadystate;
    write_xls = options_.shock_decomp.write_xls;

initial_date = options_.initial_date;

switch realtime_
    
    case 0
        z = oo_.shock_decomposition;
    
    case 1 % realtime
        if vintage_
            z = oo_.realtime_shock_decomposition.(['time_' int2str(vintage_)]);
            fig_names=[fig_names '_realtime_' char(initial_date+vintage_-1)];
        else
            z = oo_.realtime_shock_decomposition.pool;
            fig_names=[fig_names '_realtime_pool'];
        end
    
    case 2 % conditional
        if vintage_
            z = oo_.conditional_shock_decomposition.(['time_' int2str(vintage_)]);
            initial_date = options_.initial_date+vintage_-forecast_;
            fig_names=[fig_names '_conditional_' int2str(forecast_) 'step_' char(initial_date)];
        else
            z = oo_.conditional_shock_decomposition.pool;
            fig_names=[fig_names '_conditional_pool'];
        end
        
    case 3 % forecast
        if vintage_
            z = oo_.realtime_forecast_shock_decomposition.(['time_' int2str(vintage_)]);
            initial_date = options_.initial_date+vintage_-1;
            fig_names=[fig_names '_forecast_' int2str(forecast_) 'step_' char(initial_date)];
        else
            z = oo_.realtime_forecast_shock_decomposition.pool;
            fig_names=[fig_names '_forecast_1step_pool'];
        end
end

gend = size(z,3);
if options_.use_shock_groups
    shock_groups = M_.shock_groups.(options_.use_shock_groups);
    shock_ind = fieldnames(shock_groups);
    ngroups = length(shock_ind);
    fig_names=[fig_names '_group_' options_.use_shock_groups];
    shock_names = shock_ind;
    for i=1:ngroups,
       shock_names{i} = (shock_groups.(shock_ind{i}).label);
    end
    zz = zeros(endo_nbr,ngroups+2,gend);
    kcum=[];
    for i=1:ngroups
        for j = shock_groups.(shock_ind{i}).shocks
            k = find(strcmp(j,cellstr(M_.exo_names)));
            zz(:,i,:) = zz(:,i,:) + z(:,k,:);
            z(:,k,:) = 0;
            kcum = [kcum k];
        end
    end
    zothers = sum(z(:,1:nshocks,:),2);
    shock_groups.(['group' int2str(ngroups+1)]).shocks =  cellstr(M_.exo_names(find(~ismember([1:M_.exo_nbr],kcum)),:))';
    M_.shock_groups.(options_.use_shock_groups)=shock_groups;
    if any(any(zothers)),
        shock_names = [shock_names; {'Others + Initial Values'}];
    end        
    zz(:,ngroups+1,:) = sum(z(:,1:nshocks+1,:),2);
    zz(:,ngroups+2,:) = z(:,nshocks+2,:);
    z = zz;
else
    shock_names = M_.exo_names;
end

        func = @(x) colorspace('RGB->Lab',x);
        MAP = distinguishable_colors(size(z,2)-1,'w',func);
%         MAP = [MAP; MAP(end,:)];
        MAP(end,:) = [0.7 0.7 0.7];
%         MAP = [MAP; [0.7 0.7 0.7]; [0.3 0.3 0.3]];

if isempty(options_.colormap),
    options_.colormap = MAP;
end
steady_state = oo_.steady_state;
fig_mode=type;

switch type

    case '' % default

    case 'qoq' 

    case 'yoy'
        z=z(:,:,1:end-3)+z(:,:,2:end-2)+z(:,:,3:end-1)+z(:,:,4:end);
        if ~isempty(initial_date),
            initial_date = initial_date+3;
        else
            initial_date = dates('0Q4');
        end
        fig_mode = type;
        steady_state = 4*steady_state;
        
    case 'aoa'

        if isempty(initial_date),
            initial_date = dates('1Y');
            t0=4;
        else
            initial_date = dates([int2str(initial_date.time(1)) 'Y']);
            t0=4-initial_date.time(2)+1;
        end
        z=z(:,:,t0:4:end);
        fig_mode = 'AoA';

    otherwise

        error('plot_shock_decomposition:: Wrong type')

end
if steadystate
    options_.shock_decomp.steady_state=steady_state;
end
options_.shock_decomp.fig_mode=fig_mode;
options_.shock_decomp.fig_names=fig_names;
if detail_plot,
    graph_decomp_detail(z,shock_names,M_.endo_names,i_var,initial_date,M_,options_)
else
    graph_decomp(z,shock_names,M_.endo_names,i_var,initial_date,M_,options_,options_.shock_decomp);
end

if write_xls
    WriteShockDecomp2Excel(z,shock_names,M_.endo_names,i_var,initial_date,M_,options_,options_.shock_decomp);
end