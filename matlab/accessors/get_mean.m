function y0 = get_mean(varargin)
% function x = get_mean(vname1, vname2, <order>)
% returns the steady-state of a variable identified by its name
%
% INPUTS:
%   vname1, vname2, ... :  list of variable names
%   order: if integer 1 or 2, optionally last input can trigger the order
%   at which steady state is computed
%
% OUTPUTS
%   x:      steady state values
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
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

global M_ oo_ options_

if ~isempty(regexp(varargin{end},'\d','ONCE')) && isempty(regexp(varargin{end},'\D','ONCE'))
    order=eval(varargin{end});
else
    order=1;
end
if order==1
    ys_ = oo_.steady_state;
    ys_ = evaluate_steady_state(ys_,M_,options_,oo_,1);
elseif order==2
    ys_ = oo_.dr.ys;
    ys_(oo_.dr.order_var)=ys_(oo_.dr.order_var)+oo_.dr.ghs2./2;
else
    return
end
lgy_ = M_.endo_names;

mfys=nan(length(varargin),1);
for j=1:length(varargin)
    mfys(j) = find(strcmp(varargin{j},lgy_));
end

y0 = ys_(mfys);
