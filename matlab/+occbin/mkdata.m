function [data_mat]=mkdata(n_periods,dr_A,dr_B,endo_names,exo_names,wish_list,shock_pos,shock_size_vec,var_init)
% function [data_mat]=mkdata(n_periods,dr_A,dr_B,endo_names,exo_names,wish_list,shock_pos,shock_size_vec,var_init)
% Inputs:
% - n_periods       [integer]   number of simulation periods
% - dr_A            [double]    [n by n] transition matrix
% - dr_B            [double]    [n by nexo] shock response matrix
% - endo_names      [cell]      name of endogenous variables
% - exo_names       [cell]      name of exogenous variables 
% - wish_list       [cell]      name of requested variables for output
% - shock_pos       [integer]   index of shocks 
% - shock_size_vec  [double]    [shock periods by 1] vector of 
% - var_init      [double]      [n by 1] vector of initial values (incl. states)
%
% Outputs:
% - data_mat        [double]    [n_periods by n] vector of regime number indices

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

% given decision rule
neqs = size(endo_names,1);

if nargin<9
    var_init = zeros(neqs,1);
end

if nargin<8
    shock_size_vec=1;
end

if nargin<7
    error('Not enough inputs')
end

history = zeros(neqs,n_periods+1);

% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history(:,1)= var_init;

lengthshock = size(shock_size_vec,1);

err_vec = zeros(size(exo_names,1),1);

for i = 2:n_periods+1
    if i<=(lengthshock+1)
        err_vec(shock_pos) = shock_size_vec(i-1,:);
        history(:,i) = dr_A * history(:,i-1)+dr_B*err_vec;
    else
        % update endogenous variables
        history(:,i) = dr_A * history(:,i-1);
    end
end

% extract desired variables
if ~isempty(wish_list)
    n_wish=size(wish_list,1);
    wish_pos = zeros(n_wish,1);    
    for i=1:n_wish
        wish_pos(i) = strmatch(wish_list(i,:),endo_names,'exact');
    end
    data_mat = history(wish_pos,2:end)';
else
    data_mat = history(:,2:end)';
end
