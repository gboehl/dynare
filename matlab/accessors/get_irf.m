
function y0 = get_irf(exo,varargin)
% function x = get_irf(exoname, vname1, vname2, ...)
% returns IRF to individual exogenous for a list of variables and adds the
% steady state
%
% INPUTS:
%   exo: exo variable name
%   vname1, vname2, ... :  list of variable names
%
% OUTPUTS
%   x:      irf matrix [time x number of variables]
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

global M_ oo_

ys_ = [oo_.steady_state];
y0=zeros(length(oo_.irfs.([varargin{1} '_' exo]))+1,length(varargin));


[i_var,nvar] = varlist_indices(varargin,M_.endo_names);


for j=1:nvar
    %     mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=[0; oo_.irfs.([ varargin{j} '_' exo ])']+ys_(i_var(j));
end
