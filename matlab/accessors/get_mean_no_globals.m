function y0 = get_mean_no_globals(M_, oo_, options_, varargin)
% function y0 = get_mean_no_globals(M_, oo_, options_, varargin)
% returns the mean of a variable identified by its name
%
% INPUTS:
%   M_                  [structure] describing the model
%   oo_                 [structure] storing the results
%   options_            [structure] describing the options
%   vargargin           inputs containing
%                       - vname1, vname2, ... :  list of variable names
%                       - order: if integer 1 or 2, optionally last input can trigger the order
%                           at which steady state is computed
%
% OUTPUTS
%   y0:      mean values
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019-2023 Dynare Team
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

if ~isempty(regexp(varargin{end},'\d','ONCE')) && isempty(regexp(varargin{end},'\D','ONCE'))
    order=eval(varargin{end});
    nvars=length(varargin)-1;
else
    order=1;
    nvars=length(varargin);
end
if order==1
    if isfield(oo_,'dr') && isfield(oo_.dr,'ys')
        ys_=oo_.dr.ys;
    else
        ys_ = oo_.steady_state;
        ys_ = evaluate_steady_state(ys_,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,true);
    end
elseif order==2
    if isfield(oo_,'dr') && isfield(oo_.dr,'ys')
        ys_=oo_.dr.ys;
        if ~isfield(oo_.dr,'ghs2')
            error('get_mean: ghs2 needs to be present in oo_ to compute mean at order=2')
        else
            ys_(oo_.dr.order_var)=ys_(oo_.dr.order_var)+oo_.dr.ghs2./2;
        end
    else
        error('get_mean: decision rules need to be present in oo_ to compute mean') 
    end
else
    error('get_mean: order>2 not implemented')
end

mfys=NaN(nvars,1);
for j=1:nvars
    endo_index=find(strcmp(varargin{j},M_.endo_names));
    if isempty(endo_index)
        error('get_mean: unknown variables %s requested',varargin{j})
    else
    mfys(j) = find(strcmp(varargin{j},M_.endo_names));
    end
end

y0 = ys_(mfys);
