var c k x;

parameters alpha gam delta beta ;

alpha=0.5;
gam=0.5;
delta=0.02;
beta=0.05;

model;
x = x(-1)*1.02;
c + k - x^(1-alpha)*k(-1)^alpha - (1-delta)*k(-1);
c^(-gam) - (1+beta)^(-1)*(alpha*x(+1)^(1-alpha)*k^(alpha-1) + 1 - delta)*c(+1)^(-gam);
end;

verbatim;

    bgp.write(M_);
    options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000000,'MaxIterations',100000,'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-6,'StepTolerance',1e-6);
    y = 1+(rand(M_.endo_nbr,1)-.5)*.5;
    g = ones(M_.endo_nbr,1);% 1+(rand(M_.endo_nbr,1)-.5)*.1;
    [y, fval, exitflag] = fsolve(@ramsey.bgpfun, [y;g], options);
    assert(max(abs(y(M_.endo_nbr+(1:M_.orig_endo_nbr))-1.02))<1e-6)

end;
