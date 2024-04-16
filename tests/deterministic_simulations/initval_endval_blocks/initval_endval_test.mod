var Y K I L q KL_ratio Profit d;

predetermined_variables K;

varexo eta;

parameters alpha beta delta psi w Pi A r;
alpha   = 1/3;     
delta   = 0.02;   
beta    = 0.98;     
psi     = 0.0003393;     
Pi      = 1; 
A       = 4.942921;    
w       = 20.2966;
r   = 0.02;

model;

//[name = 'FOC wrt L', mcp = 'L>0']
(1-alpha)*A*((K)^alpha)*(L^(-alpha)) = w*(1-eta);

[name = 'FOC wrt K(+1)']
beta*alpha*A*(K(+1)^(alpha-1))*(L(+1)^(1-alpha)) + beta*(1-delta)*q(+1) = q ⟂ q>0;

KL_ratio = K/L;

//[name = 'FOC wrt I', mcp = 'I>0']
I = (1/psi)*(1 - (Pi/q));

//[name = 'Law of Motion for K', mcp = 'K>0']
K(+1) = (1-delta)*K + I - (psi/2)*(I^2);

//[name = 'Output Function', mcp = 'Y>0']
Y = A*(K^alpha)*(L^(1-alpha));

Profit = Y - w*L - Pi*I;

//[name = 'dividend definition', mcp = 'd>0']
d = Profit;

end;

steady_state_model;
KL_ratio = ( (w*(1-eta))/((1-alpha)*A) )^(1/alpha);
q = (beta*alpha*A*(KL_ratio^(alpha-1))) / (1-beta*(1-delta));
I = (1/psi) * (1 - (Pi/q));
K = (1/delta)*(I - (psi/2)*(I^2));
L = (1/KL_ratio) * K;
Y = A*(K^(alpha))*(L^(1-alpha));
Profit = Y - w*L - Pi*I;
d = Profit;
end;

initval;
K = 5000;
end;

endval;
eta = 0;
K = 1000;
I = 1938.246;
q = 1.07481;
L = 42.43821;
KL_ratio = 233.6644;
Y = 36641.31;
end;

perfect_foresight_setup(periods=200);

if oo_.endo_simul(strmatch('K',M_.endo_names,'exact'),1)~=5000
    error('initval does not match')
end

if oo_.endo_simul(strmatch('K',M_.endo_names,'exact'),end)~=1000
    error('endval does not match')
end
