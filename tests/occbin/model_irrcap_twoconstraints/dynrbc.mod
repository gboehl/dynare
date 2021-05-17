//Tests Occbin with 2 constraints

// variables
var a, c, i, k, lambdak;    
 
// innovations to shock processes
varexo erra;


// parameters
parameters ALPHA, DELTAK, BETA, GAMMAC, RHOA, PHI, PSI, PSINEG, INEG, IRR;

model;

# zkss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
# zcss = -DELTAK*zkss + zkss^ALPHA;
# ziss = DELTAK*zkss;
# zuss = (zcss^(1-GAMMAC)-1)/(1-GAMMAC);
# zvss = zuss/(1-BETA);

/////////////////////////////////////////////////////////////////
// 1.
-exp(c)^(-GAMMAC)*(1+2*INEG*PSI*(exp(k)/exp(k(-1))-1)/exp(k(-1)))
+ BETA*exp(c(1))^(-GAMMAC)*((1-DELTAK)-2*INEG*PSI*(exp(k(1))/exp(k)-1)*
  (-exp(k(1))/exp(k)^2)+ALPHA*exp(a(1))*exp(k)^(ALPHA-1))= 
  -lambdak+BETA*(1-DELTAK)*lambdak(1);

// 2.
exp(c)+exp(k)-(1-DELTAK)*exp(k(-1))+
INEG*PSI*(exp(k)/exp(k(-1))-1)^2=exp(a)*exp(k(-1))^(ALPHA);

// 3.
exp(i) = exp(k)-(1-DELTAK)*exp(k(-1));

// 4.
lambdak*(1-IRR) + IRR*(i - log(PHI*ziss)) = 0;

// 5. 
a = RHOA*a(-1)+erra;
end;

occbin_constraints;
name 'IRR'; bind i<PHI-1; relax lambdak<0;
name 'INEG'; bind i<-0.000001;
end;

@#include "dynrbc_common.inc"
