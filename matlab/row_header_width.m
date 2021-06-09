function w=row_header_width(M_,estim_params_,bayestopt_)
% This function computes the width of the row headers for
% the estimation results
%
% INPUTS
%   estim_params_    [structure]
%   M_               [structure]
%   bayestopt_       [structure]
%
% OUTPUTS
%   w                integer
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2018 Dynare Team
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

np = estim_params_.np;
nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;

w = 0;
if np
    w = cellofchararraymaxlength(bayestopt_.name);
end
if nvx
    w = max(w, cellofchararraymaxlength(M_.endo_names(estim_params_.var_exo(1:nvx,1))));
end
if nvn
    w = max(w, cellofchararraymaxlength(M_.endo_names(estim_params_.var_endo(1:nvn,1))));
end
if ncx
    for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        w = max(w, length(M_.exo_names{k1})+length(M_.exo_names{k2}));

    end
end
if ncn
    for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        w = max(w, length(M_.endo_names{k1})+length(M_.endo_names{k2}));
    end
end
