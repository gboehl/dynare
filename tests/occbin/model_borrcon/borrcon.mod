// --+ options: nostrict +--
//Test Occbin with 1 constraint and redundant shocks; also checks whether defaults for error_* are correct

var b ${b}$ (long_name='borrowing')
    c ${c}$ (long_name='vonsumption')
    ec ${E(c_t)}$ (long_name='expected consumption')
    lb ${\lambda}$ (long_name='Lagrange multiplier')
    y ${y}$ (long_name='Output')
    c_hat ${\hat c}$
    b_hat ${\hat c}$
    y_hat ${\hat y}$    
    ;
varexo u $u$;

parameters RHO ${\rho}$, BETA ${\beta}$, M $M$, R $R$, SIGMA ${\sigma}$, GAMMAC $\gamma_c$, relax_borrcon  ;

model;
ec = c(1);
c = y + b - R*b(-1) ;
(1-relax_borrcon)*(b - M*y) +  relax_borrcon*lb = 0;
lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC ;
log(y) = RHO*log(y(-1)) + u ;
c_hat = log(c) - log(steady_state(c));
b_hat = log(b) - log(steady_state(b));
y_hat = log(y) - log(steady_state(y));
end;

occbin_constraints;
% name 'relax_borrcon'; bind lb<-lb_ss; relax b>M*y; error_bind abs(lb+lb_ss); error_relax abs(b-M*y);
name 'relax_borrcon'; bind lb<-lb_ss; relax b>M*y;
end;
    
steady_state_model;
b=M;
c=1+M-R*M;
ec=c;
lb=(1-BETA*R)/c^GAMMAC;
y=1;
end;

R = 1.05;
BETA = 0.945;
RHO   = 0.9;
SIGMA = 0.05;
M = 1;
GAMMAC = 1;
relax_borrcon = 0;

shocks;
var u; stderr SIGMA;
end;

@#include "borrcon_common.inc"