% RBC model used in replication files of 
% Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
% Adapted by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright (C) 2020 Dynare Team

var k   $K$
    c   $C$
    a   $A$
    iv  $I$
    y   $Y$
    la  $\lambda$
    n   $N$
    rk  ${r^k}$
    w   $W$
    ;
    
predetermined_variables k;

varexo u_a ${\varepsilon^{a}}$
    ;
    
parameters DELTA    $\delta$
    BETTA           $\beta$
    B               $B$
    ETAl            $\eta_l$
    ETAc            $\eta_c$
    THETA           $\theta$
    ALFA            $\alpha$
    RHOA            $\rho_a$
    ;       

DELTA           = 0.025;
BETTA           = 0.984;
B               = 0.5;
ETAl            = 1; 
ETAc            = 2; 
THETA           = 3.48;
ALFA            = 0.667;
RHOA            = 0.979;

model;    
0 = -exp(la) +(exp(c)-B*exp(c(-1)))^(-ETAc) - BETTA*B*(exp(c(+1))-B*exp(c))^(-ETAc);
0 = -THETA*(1-exp(n))^-ETAl + exp(la)*exp(w);
0 = -exp(la) + BETTA*exp(la(+1))*(exp(rk(+1)) + (1-DELTA));
0 = -exp(a)*(1-ALFA)*exp(k)^(-ALFA)*exp(n)^(ALFA) + exp(rk);
0 = -exp(a)*ALFA*exp(k)^(1-ALFA)*exp(n)^(ALFA-1) + exp(w);
0 = -exp(c) - exp(iv) + exp(y);
0 = -exp(y) + exp(a)*exp(k)^(1-ALFA)*exp(n)^(ALFA);
0 = -exp(k(+1)) + (1-DELTA)*exp(k) + exp(iv);
0 = -log(exp(a)) + RHOA*log(exp(a(-1))) + u_a;
end;

steady_state_model;
A = 1;
RK = 1/BETTA - (1-DELTA);
K_O_N = (RK/(A*(1-ALFA)))^(-1/ALFA);
W = A*ALFA*(K_O_N)^(1-ALFA);
IV_O_N = DELTA*K_O_N;
Y_O_N = A*K_O_N^(1-ALFA);
C_O_N = Y_O_N - IV_O_N;

[N, info]=RBC_MoM_steady_helper(THETA,ETAl,ETAc,BETTA,B,C_O_N,W);
C=C_O_N*N;
Y=Y_O_N*N;
IV=IV_O_N*N;
K=K_O_N*N;
LA = (C-B*C)^(-ETAc)-BETTA*B*(C-B*C)^(-ETAc);

k=log(K);
c=log(C);
a=log(A);
iv=log(IV);
y=log(Y);
la=log(LA);
n=log(N);
rk=log(RK);
w=log(W);

end;