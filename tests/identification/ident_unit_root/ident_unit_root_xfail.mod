//Tests Identification command with ML and unit roots/diffuse filter option;
//Should not work because of observed unit root variable

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