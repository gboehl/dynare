function [zdata, T, R, CONST, ss, update_flag]=mkdatap_anticipated_2constraints_dyn(n_periods,DM,T_max,...
    binding_indicator,irfshock_pos,scalefactor_mod,init,update_flag)
% function [zdata, T, R, CONST, ss, update_flag]=mkdatap_anticipated_2constraints_dyn(n_periods,DM,T_max,...
%     binding_indicator,irfshock_pos,scalefactor_mod,init,update_flag)
%
% Inputs:
% - n_periods           [double]        number for periods for simulation
% - DM                  [structure]     Dynamic model 
% - T_max               [Tmax]          last period where constraints bind
% - binding_indicator   [T+1]           indicator for constraint violations      
% - irfshock_pos        [double]        shock position 
% - scalefactor_mod     [double]        shock values
% - init                [double]        [N by 1] initial value of endogenous variables
% - update_flag         [boolean]       flag whether to update results
%
% Output:
% - zdata               [double]        T+1 by N matrix of simulated data
% - T                   [N by N]        transition matrix of state space
% - R                   [N by N_exo]    shock impact matrix of state space
% - CONST               [N by 1]        constant of state space
% - ss                  [structure]     state space system
% - update_flag         [boolean]       flag that results have been updated
% Original authors: Luca Guerrieri and Matteo Iacoviello 
% Original file downloaded from:
% https://www.matteoiacoviello.com/research_files/occbin_20140630.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Luca Guerrieri and Matteo Iacoviello (2015): "OccBin: A toolkit for solving
% dynamic models with occasionally binding constraints easily"
% Journal of Monetary Economics 70, 22-38

persistent dictionary

if update_flag
    dictionary=[];
    update_flag=false;
end

n_vars = DM.n_vars;
T = DM.decrulea;
CONST = zeros(n_vars,1);
R = DM.decruleb;

if nargin<7 || isempty(init)
    init=zeros(n_vars,1);
end

if nargin<6
    scalefactor_mod=1;
end

n_exo=DM.n_exo;

% Tmax = max([regimestart1(nregimes1) regimestart2(nregimes2)])-1;  % Tmax is the position of the last period
% when the constraint binds

if ~isempty(dictionary) 
    if (length(binding_indicator(:))>size(dictionary.binding_indicator,1))
        nviol_old = size(dictionary.binding_indicator,1)/2;
        tmp = zeros(length(binding_indicator)-nviol_old,size(dictionary.binding_indicator,2));
        dictionary.binding_indicator = [dictionary.binding_indicator(1:nviol_old,:); tmp; dictionary.binding_indicator(1+nviol_old:2*nviol_old,:); tmp];
    end
    if (length(binding_indicator(:))<size(dictionary.binding_indicator,1))
        binding_indicator = [binding_indicator; zeros(size(dictionary.binding_indicator,1)/2-size(binding_indicator,1),2) ];
    end
end

if T_max > 0
    
    if isempty(dictionary)
        tmp = [binding_indicator(T_max,:); zeros(n_periods,2)];
        dictionary.binding_indicator(:,1) = tmp(:);
        if (binding_indicator(T_max,1) && ~binding_indicator(T_max,2))
            temp = -(DM.Abarmat10*DM.decrulea+DM.Bbarmat10)\[DM.Cbarmat10 DM.Jbarmat10 DM.Dbarmat10];
            dictionary.ss(1).T = temp(:,1:n_vars);
            dictionary.ss(1).R = temp(:,n_vars+1:n_vars+n_exo);
            dictionary.ss(1).C = temp(:,n_vars+n_exo+1:end);
        elseif (binding_indicator(T_max,1) && binding_indicator(T_max,2))
            temp = -(DM.Abarmat11*DM.decrulea+DM.Bbarmat11)\[DM.Cbarmat11 DM.Jbarmat11 DM.Dbarmat11];
            dictionary.ss(1).T = temp(:,1:n_vars);
            dictionary.ss(1).R = temp(:,n_vars+1:n_vars+n_exo);
            dictionary.ss(1).C = temp(:,n_vars+n_exo+1:end);
        else
            temp = -(DM.Abarmat01*DM.decrulea+DM.Bbarmat01)\[DM.Cbarmat01 DM.Jbarmat01 DM.Dbarmat01];
            dictionary.ss(1).T = temp(:,1:n_vars);
            dictionary.ss(1).R = temp(:,n_vars+1:n_vars+n_exo);
            dictionary.ss(1).C = temp(:,n_vars+n_exo+1:end);
        end
        % we do not know what is the last binding regime between 10 01 and 11!\
        ireg(T_max)=1;
        icount = 1;
    else
        icount=length(dictionary.ss);
        % check if last binding regime was already stored
        tmp = 0*binding_indicator;
        tmp(1:end-T_max+1,:) = binding_indicator(T_max:end,:);
        itmp = find(~any(dictionary.binding_indicator(1:length(tmp)*2,:)-tmp(:)));
        if ~isempty(itmp)
            ireg(T_max) = itmp;
        else
            icount=icount+1;
            ireg(T_max) = icount;
            tmp = [binding_indicator(T_max,:); zeros(size(binding_indicator,1)-1,2)];
            dictionary.binding_indicator(:,icount) = tmp(:);
            if (binding_indicator(T_max,1) && ~binding_indicator(T_max,2))
                temp = -(DM.Abarmat10*DM.decrulea+DM.Bbarmat10)\[DM.Cbarmat10 DM.Jbarmat10 DM.Dbarmat10];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            elseif (binding_indicator(T_max,1) && binding_indicator(T_max,2))
                temp = -(DM.Abarmat11*DM.decrulea+DM.Bbarmat11)\[DM.Cbarmat11 DM.Jbarmat11 DM.Dbarmat11];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            else
                temp = -(DM.Abarmat01*DM.decrulea+DM.Bbarmat01)\[DM.Cbarmat01 DM.Jbarmat01 DM.Dbarmat01];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            end
        end

    end
    
    
    for i = T_max-1:-1:1        
        tmp = 0*binding_indicator;
        tmp(1:end-i+1,:) = binding_indicator(i:end,:);
        itmp = find(~any(dictionary.binding_indicator(1:length(tmp)*2,:)-tmp(:)));
        if ~isempty(itmp)
            ireg(i) = itmp;
        else
            icount=icount+1;
            ireg(i) = icount;
            dictionary.binding_indicator(1:length(tmp)*2,icount) = tmp(:);
            if (binding_indicator(i,1) && ~binding_indicator(i,2))
                temp = -(DM.Bbarmat10+DM.Abarmat10*dictionary.ss(ireg(i+1)).T)\[DM.Cbarmat10 DM.Jbarmat10 DM.Abarmat10*dictionary.ss(ireg(i+1)).C+DM.Dbarmat10];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            elseif (~binding_indicator(i,1) && binding_indicator(i,2))
                temp = -(DM.Bbarmat01+DM.Abarmat01*dictionary.ss(ireg(i+1)).T)\[DM.Cbarmat01 DM.Jbarmat01 DM.Abarmat01*dictionary.ss(ireg(i+1)).C+DM.Dbarmat01];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            elseif (binding_indicator(i,1) && binding_indicator(i,2))
                temp = -(DM.Bbarmat11+DM.Abarmat11*dictionary.ss(ireg(i+1)).T)\[DM.Cbarmat11 DM.Jbarmat11 DM.Abarmat11*dictionary.ss(ireg(i+1)).C+DM.Dbarmat11];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            else
                temp = -(DM.Bbarmat+DM.Abarmat*dictionary.ss(ireg(i+1)).T)\[DM.Cbarmat DM.Jbarmat DM.Abarmat*dictionary.ss(ireg(i+1)).C];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            end
        end
    end

    E = dictionary.ss(ireg(1)).R;
    ss = dictionary.ss(ireg(1:T_max));
else
    ss = [];    
end

% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history = zeros(n_vars,n_periods+1);
history(:,1) = init;
errvec = zeros(n_exo,1);

errvec(irfshock_pos) = scalefactor_mod;

% deal with shocks
irfpos =1;
if irfpos <=T_max
    history(:,irfpos+1) = dictionary.ss(ireg(irfpos)).T* history(:,irfpos)+...
        dictionary.ss(ireg(irfpos)).C + E*errvec;
    T = dictionary.ss(ireg(irfpos)).T;
    CONST = dictionary.ss(ireg(irfpos)).C;
    R = E;
else
    history(:,irfpos+1) = DM.decrulea*history(:,irfpos)+DM.decruleb*errvec;
end

% all other periods
for irfpos=2:n_periods+1
    if irfpos <=T_max
        history(:,irfpos+1) = dictionary.ss(ireg(irfpos)).T* history(:,irfpos)+...
            dictionary.ss(ireg(irfpos)).C;
    else
        history(:,irfpos+1) = DM.decrulea*history(:,irfpos);
    end
end

zdata = history(:,2:end)';
