function [zdata, T, R, CONST, ss, update_flag]=mkdatap_anticipated_dyn(n_periods,DM,...
    T_max,binding_indicator,irfshock_pos,scalefactor_mod,init,update_flag)
% function [zdata, T, R, CONST, ss]=mkdatap_anticipated_dyn(nperiods,DM,...
%     Tmax,binding_indicator,irfshockpos,scalefactormod,init,update_flag)
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
%
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

%Initialize outputs
n_vars = DM.n_vars;
n_exo = DM.n_exo;
T = DM.decrulea;
CONST = zeros(n_vars,1);
R = DM.decruleb;

if nargin<7 || isempty(init)
    init=zeros(n_vars,1);
end

if nargin<6
    scalefactor_mod=1;
end

 
% % get the time-dependent decision rules

if ~isempty(dictionary) 
    if (length(binding_indicator)>size(dictionary.binding_indicator,1))
        dictionary.binding_indicator = [dictionary.binding_indicator; zeros(length(binding_indicator)-size(dictionary.binding_indicator,1),size(dictionary.binding_indicator,2))];
    end
    if (length(binding_indicator(:))<size(dictionary.binding_indicator,1))
        binding_indicator = [binding_indicator; zeros(size(dictionary.binding_indicator,1)-size(binding_indicator,1),1) ];
    end
end
        
if T_max > 0    
    if isempty(dictionary)        
        temp = -(DM.Astarbarmat*DM.decrulea+DM.Bstarbarmat)\[DM.Cstarbarmat DM.Jstarbarmat DM.Dstarbarmat];
        dictionary.binding_indicator(:,1) = [1; zeros(n_periods,1)];
        dictionary.ss(1).T = temp(:,1:n_vars);
        dictionary.ss(1).R = temp(:,n_vars+1:n_vars+n_exo);
        dictionary.ss(1).C = temp(:,n_vars+n_exo+1:end);
    end
    ireg(T_max)=1;
    
    % equivalent to pre-multiplying by the inverse above if the target
    % matrix is invertible. Otherwise it yields the minimum state solution
    %P(:,:,Tmax) = -(Astarbarmat*decrulea+Bstarbarmat)\Cstarbarmat;
    %D(:,Tmax) = -(Astarbarmat*decrulea+Bstarbarmat)\Dstarbarmat;
    
    icount=length(dictionary.ss);
    
    for i = T_max-1:-1:1
        
        tmp = 0*binding_indicator;
        tmp(1:end-i+1) = binding_indicator(i:end);
        itmp = find(~any(dictionary.binding_indicator-tmp));
        if ~isempty(itmp)
            ireg(i) = itmp;
        else
            icount=icount+1;
            ireg(i) = icount;
            dictionary.binding_indicator(1:length(tmp),icount) = tmp;
            if binding_indicator(i)
                temp = -(DM.Bstarbarmat+DM.Astarbarmat*dictionary.ss(ireg(i+1)).T)\[DM.Cstarbarmat DM.Jstarbarmat DM.Astarbarmat*dictionary.ss(ireg(i+1)).C+DM.Dstarbarmat];
                dictionary.ss(icount).T = temp(:,1:n_vars);
                dictionary.ss(icount).R = temp(:,n_vars+1:n_vars+n_exo);
                dictionary.ss(icount).C = temp(:,n_vars+n_exo+1:end);
            else
                temp = -(DM.Bbarmat+DM.Abarmat*dictionary.ss(ireg(i+1)).T)\[DM.Cbarmat DM.Jbarmat (DM.Abarmat*dictionary.ss(ireg(i+1)).C)];
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

% deal with predetermined conditions
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
