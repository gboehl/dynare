% check for correct error message if token cannot be interpreted
// variables
var a, c, i, k, lambdak;    
 
// innovations to shock processes
varexo erra;


// parameters
parameters ALPHA, DELTAK, BETA, GAMMAC, RHOA, PHI, PSI, PSINEG, INEG, IRR;

model(occbin);

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
[pswitch = 'INEG', 
//         bind = 'exp(i+i_ss)<-0.000001', 
//         relax = 'exp(i+i_ss)>-0.000001' ]
        bind = 'i<-b', 
        relax = 'i>-0.000001' ]
exp(i) = exp(k)-(1-DELTAK)*exp(k(-1));

// 4.
[pswitch = 'IRR', 
        bind = 'i<PHI-1', 
        relax = 'lambdak<0' ]
lambdak*(1-IRR) + IRR*(i - log(PHI*ziss)) = 0;

// 5. 
a = RHOA*a(-1)+erra;


end;

steady_state_model;
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;


k = log(kss);
c = log(css);
i = log(iss);
lambdak = 0;
a=0;
end;

BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;
GAMMAC=2;
RHOA = 0.9;
PHI = 0.975;
PSI = 5;        % adjustment cost for capital if investment is negative
INEG = 0;
IRR = 0;

shocks;
  var erra; stderr 0.015;
end;

steady;

stoch_simul(order=1,nocorr,nomoments,irf=0,noprint);
