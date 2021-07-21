//Tests Occbin with 2 constraints

// variables
var a, c, i, k, lambdak;    
 
// innovations to shock processes
varexo erra;


// parameters
parameters ALPHA, DELTAK, BETA, GAMMAC, RHOA, PHI, PSI;

model;

# zkss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
# zcss = -DELTAK*zkss + zkss^ALPHA;
# ziss = DELTAK*zkss;
# zuss = (zcss^(1-GAMMAC)-1)/(1-GAMMAC);
# zvss = zuss/(1-BETA);

/////////////////////////////////////////////////////////////////
// 1.
[name='Euler', bind = 'INEG'] 
-exp(c)^(-GAMMAC)*(1+2*PSI*(exp(k)/exp(k(-1))-1)/exp(k(-1)))
+ BETA*exp(c(1))^(-GAMMAC)*((1-DELTAK)-2*PSI*(exp(k(1))/exp(k)-1)*
  (-exp(k(1))/exp(k)^2)+ALPHA*exp(a(1))*exp(k)^(ALPHA-1))= 
  -lambdak+BETA*(1-DELTAK)*lambdak(1);

[name='Euler', relax = 'INEG']
-exp(c)^(-GAMMAC) + BETA*exp(c(1))^(-GAMMAC)*(1-DELTAK+ALPHA*exp(a(1))*exp(k)^(ALPHA-1))= 
  -lambdak+BETA*(1-DELTAK)*lambdak(1);

// 2.
[name='Budget constraint',bind = 'INEG'] 
exp(c)+exp(k)-(1-DELTAK)*exp(k(-1))+PSI*(exp(k)/exp(k(-1))-1)^2=exp(a)*exp(k(-1))^(ALPHA);

[name='Budget constraint',relax = 'INEG']
exp(c)+exp(k)-(1-DELTAK)*exp(k(-1))=exp(a)*exp(k(-1))^(ALPHA);

// 3.
exp(i) = exp(k)-(1-DELTAK)*exp(k(-1));

// 4.
[name='investment',bind='IRR,INEG']
(i - log(PHI*ziss)) = 0;
[name='investment',relax='IRR']
lambdak=0;
[name='investment',bind='IRR',relax='INEG']
(i - log(PHI*ziss)) = 0;

// 5. 
a = RHOA*a(-1)+erra;
end;

occbin_constraints;
name 'IRR'; bind i-steady_state(i)<log(PHI); relax lambdak<0;
name 'INEG'; bind i-steady_state(i)<-0.000001;
end;

@#include "dynrbc_common.inc"
