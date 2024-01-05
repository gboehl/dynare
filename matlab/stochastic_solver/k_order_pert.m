function [dr,info] = k_order_pert(dr,M_,options_)
% Compute decision rules using the k-order DLL from Dynare++

% Copyright © 2009-2024 Dynare Team
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

info = 0;

order = options_.order;

if order>1 && options_.loglinear
    error('The loglinear-option currently only works at order 1')
end
if M_.maximum_endo_lead == 0 && order>1
    error(['2nd and 3rd order approximation not implemented for purely ' ...
           'backward models'])
end

if options_.aim_solver
    error('Option aim_solver is not compatible with k_order_solver')
end
if options_.dr_cycle_reduction
    error('Option dr=cycle_reduction is not compatible with k_order_solver')
end
if options_.dr_logarithmic_reduction
    error('Option dr=logarithmic_reduction is not compatible with k_order_solver')
end

try
    [dynpp_derivs, dyn_derivs] = k_order_perturbation(dr,M_,options_);
catch ME
    disp(ME.message)
    info(1)=9;
    return
end

% Copy Dynare++ tensors

for i = 0:order
    gname = [ 'g_' num2str(i) ];
    dr.(gname) = dynpp_derivs.(gname);
end
if options_.pruning
   dr.pruning = dynpp_derivs.pruning;
end


% Fill equivalent Dynare matrices (up to 3rd order)

dr.ghx = dyn_derivs.gy;
dr.ghu = dyn_derivs.gu;

if options_.loglinear
    k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;
    dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
             repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
    dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
end

if order >= 2
    dr.ghxx = dyn_derivs.gyy;
    dr.ghxu = dyn_derivs.gyu;
    dr.ghuu = dyn_derivs.guu;
    dr.ghs2 = dyn_derivs.gss;
end

if order >= 3
    dr.ghxxx = dyn_derivs.gyyy;
    dr.ghxxu = dyn_derivs.gyyu;
    dr.ghxuu = dyn_derivs.gyuu;
    dr.ghuuu = dyn_derivs.guuu;
    dr.ghxss = dyn_derivs.gyss;
    dr.ghuss = dyn_derivs.guss;
end
