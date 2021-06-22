function [dr, info, M_, oo_]=discretionary_policy_1(Instruments, M_, options_, oo_)
% Higher-level function for solving discretionary optimal policy
% INPUTS
% - Instruments   [cell]          array containing instrument names
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% - oo_           [structure]     Matlab's structure containing the results (oo_).
%
% OUTPUTS
% - dr            [structure]     Reduced form model.
% - info          [integer]       scalar or vector, error code.
% - M_            [structure]     Matlab's structure describing the model (M_).
% - oo_           [structure]     Matlab's structure containing the results (oo_).

% Copyright (C) 2007-2020 Dynare Team
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

persistent Hold

info = 0;

dr=oo_.dr; %initialize output argument

beta = get_optimal_policy_discount_factor(M_.params, M_.param_names);

%call steady_state_file if present to update parameters
if options_.steadystate_flag
    % explicit steady state file
    [ys,M_.params,info] = evaluate_steady_state_file(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_, ...
                                                    options_,false);
    if info(1)
        return;
    end
else
    ys=zeros(M_.endo_nbr,1);
end
[U,Uy,W] = feval([M_.fname,'.objective.static'],zeros(M_.endo_nbr,1),[], M_.params);
if any(any(isnan(Uy)))
    info = 64 ; %the derivatives of the objective function contain NaN
    return;
end
if any(any(Uy~=0))
    if options_.debug
        non_zero_derivs=find(any(Uy~=0));
        for ii=1:length(non_zero_derivs)
            non_zero_deriv_names{ii,1} = M_.endo_names{non_zero_derivs(ii)};
        end
        disp_string=[non_zero_deriv_names{1,:}];
        for ii=2:size(non_zero_deriv_names,1)
            disp_string=[disp_string,', ',non_zero_deriv_names{ii,:}];
        end
        fprintf('\nThe derivative of the objective function w.r.t. to variable(s) %s is not 0\n',disp_string);
    end
    info = 66;
    return;
end

W=reshape(W,M_.endo_nbr,M_.endo_nbr);

klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv=M_.lead_lag_incidence';
% Find the jacobian
z = repmat(ys,1,klen);
iyr0 = find(iyv(:)) ;

z = z(iyr0);
it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

[junk,jacobia_] = feval([M_.fname '.dynamic'],z, [zeros(size(oo_.exo_simul)) ...
                    oo_.exo_det_simul], M_.params, ys, it_);
if max(abs(junk))>options_.solve_tolf
     info = 65; %the model must be written in deviation form and not have constant terms or have a steady state provided
     return;
end

Indices={'lag','contemp','lead'};
iter=1;
for j=1:numel(Indices)
    A.(Indices{j})=zeros(M_.eq_nbr,M_.endo_nbr);
    if strcmp(Indices{j},'contemp')||(strcmp(Indices{j},'lag') && M_.maximum_lag)||(strcmp(Indices{j},'lead') && M_.maximum_lead)
        [~,row,col]=find(M_.lead_lag_incidence(iter,:));
        A.(Indices{j})(:,row)=jacobia_(:,col);
        iter=iter+1;
    end
end
B=jacobia_(:,nnz(iyv)+1:end);

%%% MAIN ENGINE %%%

if ~isempty(Hold)
    [H,G,info]=discretionary_policy_engine(A.lag,A.contemp,A.lead,B,W,M_.instr_id,beta,options_.dp.maxit,options_.discretionary_tol,options_.qz_criterium,Hold);
else
    [H,G,info]=discretionary_policy_engine(A.lag,A.contemp,A.lead,B,W,M_.instr_id,beta,options_.dp.maxit,options_.discretionary_tol,options_.qz_criterium);
end

if info
    return
else
    Hold=H; %save previous solution
            % Hold=[]; use this line if persistent command is not used.
end

%write back solution to dr
dr.ys =ys;
dr=set_state_space(dr,M_,options_);
T=H(dr.order_var,dr.order_var);
dr.ghu=G(dr.order_var,:);
if M_.maximum_endo_lag
    Selection=M_.lead_lag_incidence(1,dr.order_var)>0;%select state variables
end
dr.ghx=T(:,Selection);
oo_.dr = dr;
