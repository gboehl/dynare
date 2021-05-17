function oo_ = shock_decomposition(oo_, M_, options_, vname)
% function oo_ = shock_decomposition(oo_, M_, options_, vname)
%
% INPUTS
% - oo_           [structure]     Matlab's structure containing the results (oo_).
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% - vname         [cell]          array of variable names
% OUTPUT
% - oo_           [structure]     Matlab's structure containing the results (oo_).

% Copyright (C) 2021 Dynare Team
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


% list of options 
% T0, use_shock_groups, vname, file_name, nfrcst, init_names_

shock_decomp_options = options_.occbin.shock_decomp;

if nargin<4 || isempty(vname)
    vname=cellstr(M_.endo_names(1:M_.orig_endo_nbr,:));
end

opts_simul = options_.occbin.simul;
opts_simul.periods = shock_decomp_options.periods;
opts_simul.check_ahead_periods = shock_decomp_options.check_ahead_periods;
opts_simul.maxit = shock_decomp_options.maxit;
use_shock_groups = shock_decomp_options.use_shock_groups;

init_names_ = shock_decomp_options.init_names_;
nfrcst = shock_decomp_options.nfrcst;

%% new dynare grouping
if isempty(use_shock_groups)
    use_shock_groups = 'ALL';
    ngroups = M_.exo_nbr;
    ex_names_ = cell(ngroups,1);
    for i=1:ngroups
        ex_names_{i} = M_.exo_names(i);
    end
    shock_decomp_options.main_effect=0;
    shock_decomp_options.additive = 0;
else
    shock_groups = M_.shock_groups.(use_shock_groups);
    shock_ind = fieldnames(shock_groups);
    ngroups = length(shock_ind);
    ex_names_ = shock_ind;
    for i=1:ngroups
        ex_names_{i} = shock_groups.(shock_ind{i}).shocks;
    end
end
%%  

if iscell(vname{1})
    vname0=vname;
    clear vname
    vname = cell(1,length(vname0));
    decomp_type = cell(1,length(vname0));
    gtrend = cell(1,length(vname0));
    var_type = nan(1,length(vname0));
    for j=1:length(vname0)
        vname{j}=vname0{j}{1};
        decomp_type{j}=vname0{j}{2};
        if strcmpi(decomp_type{j},'aoa')
            gtrend{j}=vname0{j}{3};
            if strcmpi(vname0{j}{4},'flow')
            var_type(j)=1;
            elseif strcmpi(vname0{j}{4},'deflator')
            var_type(j)=2;
            elseif strcmpi(vname0{j}{4},'stock')
            var_type(j)=0;
            else
                error('wrong var type for aoa decomp')
            end
        end
    end
    clear vname0
else
    decomp_type = cell(1,length(vname));
    for j=1:length(vname)
        decomp_type{j}='qoq';
    end
end

% set time-varying state space matrices
TT = oo_.occbin.smoother.T0;
RR = oo_.occbin.smoother.R0;
CC = oo_.occbin.smoother.C0;

load (options_.datafile, 'T');

indx_init=zeros(length(init_names_),1);
for j=1:length(init_names_)
    indx_init(j)=strmatch(init_names_{j},M_.endo_names(oo_.dr.order_var,:),'exact');
end

if shock_decomp_options.debug
    % re-compute smoothed variables from T, R and  smoothed shocks etahat
    a1=oo_.occbin.smoother.alphahat(:,1);
    as=a1;
    as0=zeros(length(a1),1);
    as0(indx_init)=a1(indx_init);
    etahat=[oo_.occbin.smoother.etahat zeros(size(oo_.occbin.smoother.etahat,1),nfrcst)];
    for j=1:size(etahat,2)-1
        TM = TT(:,:,j+1);
        RM = RR(:,:,j+1);
        CONST = CC(:,j+1);
        as(:,j+1)=TM*as(:,j)+RM*etahat(:,j+1)+CONST;
        as0(:,j+1)=TM*as0(:,j)+CONST;
    end
    err = max(max(abs(oo_.occbin.smoother.alphahat-as)));
    if err>1e-8
        disp('WARNING CHECK SMOOTHER:')
        disp(['simulated model with smoothed shocks differs from stored smoother by ' num2str(err)])
    else
        disp('CHECK SMOOTHER OK:')
        disp(['simulated model with smoothed shocks differs from stored smoother by less than 1e-8 (err=' num2str(err) ')'])
    end
else
    etahat=[oo_.occbin.smoother.etahat zeros(size(oo_.occbin.smoother.etahat,1),nfrcst)];
    as = oo_.occbin.smoother.alphahat;
end
gend=size(etahat,2);

%%
TT= 0:0.25:ceil(gend/4+1);
TT=TT(1:gend);
TT1= dates('0Q1'):(dates('0Q1')+gend-1);
if exist('T','var')
    TT=T(options_.first_obs)+TT;
    TT1=TT1+T(options_.first_obs)*4;
end

shock_decomp_options = set_default_option(shock_decomp_options,'TINIT',TT1(1));
tinit = max([1,find(TT1==shock_decomp_options.TINIT)]);
TINIT = char(TT1(tinit));
shock_decomp_options = set_default_option(shock_decomp_options,'file_name',['_' use_shock_groups '_' TINIT]);
file_name = shock_decomp_options.file_name;

%%

%%%%%%%%%%%%%%%%%%%% LINEAR shock decomp
as_lin=zeros(M_.endo_nbr,length(tinit:gend)); % linear smoother reconstructed
att=zeros(M_.endo_nbr,length(tinit:gend)); % linear initial condition effect
inn=zeros(M_.endo_nbr,length(tinit:gend)); % linear aggreage shocks effect without att, i.e. as_lin = inn+att;
deco=zeros(M_.endo_nbr,M_.exo_nbr,length(tinit:gend)); % full decomposition into individual shocks

att(:,1)=oo_.occbin.linear_smoother.alphahat(:,tinit);
as_lin(:,1)=oo_.occbin.linear_smoother.alphahat(:,tinit);
TM = oo_.occbin.linear_smoother.T0;
RM= oo_.occbin.linear_smoother.R0;
for j=2:length(tinit:gend)
    as_lin(:,j) = TM*as_lin(:,j-1)+RM*oo_.occbin.linear_smoother.etahat(:,j+tinit-1);
    att(:,j) = TM*att(:,j-1);
    inn(:,j) = RM*oo_.occbin.linear_smoother.etahat(:,j+tinit-1);
    if j>1
        inn(:,j) = inn(:,j) +  TM*inn(:,j-1); 
    end
    for iexo=1:M_.exo_nbr
        deco(:,iexo,j) = RM(:,iexo)*oo_.occbin.linear_smoother.etahat(iexo,j+tinit-1);
        if j>1
            deco(:,iexo,j) = deco(:,iexo,j) +  TM*deco(:,iexo,j-1);
        end
        
    end
end
as_lin=as_lin(oo_.dr.inv_order_var,:); % linear smoother reconstructed
att=att(oo_.dr.inv_order_var,:); % linear initial condition effect
inn=inn(oo_.dr.inv_order_var,:); % linear aggreage shocks effect without att
deco=deco(oo_.dr.inv_order_var,:,:); % full decomposition into individual shocks
deco(:,M_.exo_nbr+1,:)=att;
deco(:,M_.exo_nbr+2,:)=as_lin;
oo_.occbin.linear_smoother.decomp=deco;
%%%%%%%%%%%%%%%%%%%% END LINEAR shock decomp

%%

%%%%%%%%%%%%%%%%%%%% piecewise COONDITIONAL shock decomp
as_p=zeros(M_.endo_nbr,length(tinit:gend)); % smoother reconstructed
att_p=zeros(M_.endo_nbr,length(tinit:gend)); % initial condition effect
inn_p=zeros(M_.endo_nbr,length(tinit:gend)); % aggreage shocks effect without att, i.e. as_lin = inn+att;
deco_p=zeros(M_.endo_nbr,M_.exo_nbr,length(tinit:gend)); % full decomposition into individual shocks
reg_p=zeros(M_.endo_nbr,length(tinit:gend)); % pure regime effect (CONST 'shocks')

att_p(:,1)=oo_.occbin.smoother.alphahat(:,tinit);
as_p(:,1)=oo_.occbin.smoother.alphahat(:,tinit);
TT = oo_.occbin.smoother.T0;
RR= oo_.occbin.smoother.R0;
CC= oo_.occbin.smoother.C0;

for j=2:length(tinit:gend)
    TM = TT(:,:,j+tinit-1);
    RM = RR(:,:,j+tinit-1);
    CONST = CC(:,j+tinit-1);
    as_p(:,j) = TM*as_p(:,j-1)+RM*oo_.occbin.smoother.etahat(:,j+tinit-1)+CONST;
    att_p(:,j) = TM*att_p(:,j-1);
    inn_p(:,j) = RM*oo_.occbin.smoother.etahat(:,j+tinit-1);
    reg_p(:,j) = TM*reg_p(:,j-1)+CONST;
    if j>1
        inn_p(:,j) = inn_p(:,j) +  TM*inn_p(:,j-1) ; 
    end
    for iexo=1:M_.exo_nbr
        deco_p(:,iexo,j) = RM(:,iexo)*oo_.occbin.smoother.etahat(iexo,j+tinit-1);
        if j>1
            deco_p(:,iexo,j) = deco_p(:,iexo,j) +  TM*deco_p(:,iexo,j-1);
        end
        
    end
end
as_p=as_p(oo_.dr.inv_order_var,:); % occbin smoother reconstructed
att_p=att_p(oo_.dr.inv_order_var,:); % occbin initial condition effect
inn_p=inn_p(oo_.dr.inv_order_var,:); % occbin aggregage shocks effect without att
deco_p=deco_p(oo_.dr.inv_order_var,:,:); % occbin full decomposition into individual shocks
reg_p=reg_p(oo_.dr.inv_order_var,:); % occbin pure regime effect (CONST 'shocks')
i_reg=strmatch('EPS_REGIME',M_.exo_names,'exact');
if ~isempty(i_reg)
deco_p(:,i_reg,:)=reg_p;
end
deco_p(:,M_.exo_nbr+1,:)=att_p;
deco_p(:,M_.exo_nbr+2,:)=as_p;
%deco_p(:,M_.exo_nbr+3,:)=as_p-reg_p;
oo_.occbin.smoother.decomp=deco_p;
rr=abs(deco_p(:,1:end-1,:));
if isequal(use_shock_groups,'ALL') && ~isempty(i_reg)
    rr(:,i_reg,:)=0;
end
ww=zeros(size(rr));
for k=1:size(rr,3)
    tmp=sum(rr(:,:,k),2);
    tmp(tmp<1.e-10)=1;
    for g=1:size(rr,2)
        ww(:,g,k) = rr(:,g,k)./tmp;
    end
end
wdeco_p=deco_p;
if ~isempty(i_reg)
    wdeco_p(:,i_reg,:)=0;
end
for k=1:size(rr,3)
    if any(any(reg_p(:,k)))
        for g=1:size(rr,2)
            wdeco_p(:,g,k) = wdeco_p(:,g,k)+reg_p(:,k).*ww(:,g,k);
        end
    end
end
oo_.occbin.smoother.wdecomp=wdeco_p;
%%%%%%%%%%%%%%%%%%%% END CONDITIONAL shock decomp
%% add here other fetures when ready
% if shock_decomp_options.conditional_only ==1
%     return
% end
