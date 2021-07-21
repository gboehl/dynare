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

parameters RHO ${\rho}$, BETA ${\beta}$, M $M$, R $R$, SIGMA ${\sigma}$, GAMMAC $\gamma_c$;

model;
ec = c(1);
c = y + b - R*b(-1) ;
[name = 'borrowing', bind='borrcon']
lb = 0;
[name = 'borrowing', relax='borrcon']
b = M*y;  
lb = 1/c^GAMMAC - BETA*R/c(+1)^GAMMAC ;
log(y) = RHO*log(y(-1)) + u ;
c_hat = log(c) - log(steady_state(c));
b_hat = log(b) - log(steady_state(b));
y_hat = log(y) - log(steady_state(y));
end;

occbin_constraints;
 name 'borrcon'; bind lb<-STEADY_STATE(lb); relax b>M*y; error_bind abs(lb+STEADY_STATE(lb)); error_relax abs(b-M*y);
%name 'borrcon'; bind lb<-STEADY_STATE(lb); relax b>M*y;
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

shocks;
var u; stderr SIGMA;
end;

@#include "borrcon_common.inc"
