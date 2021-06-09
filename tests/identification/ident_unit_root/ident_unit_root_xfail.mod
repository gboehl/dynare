% Tests Identification command with ML and unit roots/diffuse filter option;
% Should not work because of observed unit root variable
% Created by Johannes Pfeifer (@JohannesPfeifer, jpfeifer@gmx.de)
% =========================================================================
% Copyright (C) 2015-2020 Dynare Team
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
% =========================================================================
var y x z delta_y;

varexo eps_x eps_z;

parameters rho sigma_z sigma_x;

// set parameter values
sigma_z=0.001;
sigma_x=0.01;
rho=0.9;

model;
z=rho*z(-1)+sigma_z*eps_z;
x=x(-1)+sigma_x*eps_x;
y=x+z;
delta_y=y-y(-1);
end;

steady_state_model;
x=0;
z=0;
y=0;
delta_y=0;
end;

//set shock variances
shocks;
    var eps_z=1;
    var eps_x=1;
end;

steady;
check;
varobs y delta_y; 
stoch_simul(order=1,irf=0);


estimated_params;
rho, 0.9;
sigma_z, 0.01;
sigma_x, 0.01;
end;
identification(diffuse_filter,advanced=1,prior_trunc=0);