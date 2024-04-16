var y pi rs rl rs_shadow rs_1;
varexo e_y e_pi;

parameters alpha_1 alpha_2 alpha_3 beta_1 beta_2 gamma_r gamma_pi gamma_y pi_bar rr_bar;


alpha_1 = 0.7;
alpha_2 = 0.99;
alpha_3 = 0.1;
beta_1 = 0.7;
beta_2 = 0.08;
gamma_r = 0.3;
gamma_pi = 2;
gamma_y = 1;
pi_bar = 2.0;
rr_bar = 1.0;

model(linear,use_dll);
    pi = pi_bar*(1-alpha_1-(1-alpha_1)*alpha_2) + alpha_1*pi(-1) + (1-alpha_1)*alpha_2*pi(+1) + alpha_3*y + e_pi;
    y = beta_1*y(-1) + (1-beta_1)*y(+1) - beta_2*(rl-rr_bar - pi(+1)) + e_y;
    rs_shadow = gamma_r*rs(-1) + (1-gamma_r)*(rr_bar + pi_bar + gamma_pi*(pi-pi_bar) + gamma_y*y);
    rs_1 = rs_shadow ⟂ rs_shadow > 0;
//    [mcp = 'pi > 1.4']
    rs = rs_1;
    rl = 400*(((1+rs/400)
    @#for i in [1:39]
        *(1+rs(+@{i})/400)
    @#endfor
    )^(1/40) - 1);
end;

steady_state_model;
y = 0;
pi = pi_bar;
rs = rr_bar+pi_bar;
rs_shadow = rs;
rs_1 = rs;
rl = rs;
end;

steady;

check;

shocks;
var e_y;
periods 1;
//values -1.226;
values -1;
end;

options_.solve_algo = 12;
options_.mcp = 1;
//options_.linear_approximation = 1;
perfect_foresight_setup(periods=200);    

perfect_foresight_solver(stack_solve_algo=7);

rplot rs;
//rplot rl;
//rplot y;
//rplot pi;
