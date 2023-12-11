function y0 = get_irf(exoname,varargin)
% function x = get_irf(exoname, varargin)
% returns IRF to individual exogenous for a list of variables and adds the
% steady state
%
% INPUTS:
%   exoname: exo variable name
%   vname1, vname2, ... :  list of variable names
%
% OUTPUTS
%   y0:      irf matrix [time x number of variables]
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

global M_ oo_

if isfield(oo_,'irfs')
    irf_fields=fieldnames(oo_.irfs);
else
    error('get_irf: No IRFs detected in oo_')
end

exo_matches=find(endsWith(irf_fields,['_' exoname]));
if isempty(exo_matches)
    error('get_irf: No IRFs for shock %s detected in oo_',exoname)
else
    if nargin>1
        endo_cell={};
        i_var=[];
        for var_iter=1:length(varargin)
            temp=startsWith(irf_fields(exo_matches),varargin{var_iter});
            if isempty(temp)
                fprintf('get_irf: No IRF for variable %s detected in oo_',varargin{var_iter})
            else
                endo_cell=[endo_cell,varargin{var_iter}];
                i_var=[i_var,strmatch(varargin(var_iter),M_.endo_names,'exact')];
            end
        end
    else
        endo_cell={};
        i_var=[];
        for var_iter=1:length(exo_matches)
            endo_cell=[endo_cell,irf_fields{var_iter}(1:end-length(exoname)-1)];
            i_var=[i_var,strmatch(endo_cell(end),M_.endo_names,'exact')];
        end
    end
end

ys_ = [oo_.steady_state];
nvars=length(endo_cell);
y0=zeros(length(oo_.irfs.([ endo_cell{1} '_' exoname ]))+1,nvars);

for j=1:nvars
    y0(:,j)=[0; oo_.irfs.([ endo_cell{j} '_' exoname ])']+ys_(i_var(j));
end


