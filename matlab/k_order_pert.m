function [dr,info] = k_order_pert(dr,M,options)
% Compute decision rules using the k-order DLL from Dynare++

% Copyright (C) 2009-2020 Dynare Team
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

order = options.order;

if order>1 && options.loglinear
    error('The loglinear-option currently only works at order 1')
end
if M.maximum_endo_lead == 0 && order>1
    error(['2nd and 3rd order approximation not implemented for purely ' ...
           'backward models'])
end

try
    [dynpp_derivs, dyn_derivs] = k_order_perturbation(dr,M,options);
catch
    info(1)=9;
    return
end

% Copy Dynare++ tensors

for i = 0:order
    gname = [ 'g_' num2str(i) ];
    dr.(gname) = dynpp_derivs.(gname);
end

% Fill equivalent Dynare matrices (up to 3rd order)

dr.ghx = dyn_derivs.gy;
dr.ghu = dyn_derivs.gu;

if options.loglinear
    k = find(dr.kstate(:,2) <= M.maximum_endo_lag+1);
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
