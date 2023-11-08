function dy = back_subst_lbj(c,ny,iyf,periods)
% Computes backward substitution in LBJ
%
% INPUTS
%    c:              matrix containing the D and d for Sébastien’s presentation
%    ny:             number of endogenous variables
%    iyf:            indices of forward variables inside the list of all endogenous variables
%    periods:        number of simulation periods
%
% OUTPUTS
%    dy:             vector of backsubstitution results

% Copyright © 2023 Dynare Team
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

ir = (periods-2)*ny+(1:ny);
irf = iyf+(periods-1)*ny;
icf = 1:size(iyf,2);
nrc = length(iyf)+1;

for i = 2:periods
    c(ir,nrc) = c(ir,nrc)-c(ir,icf)*c(irf,nrc);
    ir = ir-ny;
    irf = irf-ny;
end

dy = reshape(c(:,nrc), ny, periods);

end
