var y i pi ;

parameters a1 a2 a3 a4 a5;

a1 = -.5;
a2 =  .1;
a3 =  .9;
a4 = 1.5;
a5 = 0.5;

model;

y = y(1)*(i/pi(1))^a1;

pi = (y^a2)*(pi(1)^a3);

i = (pi^a4)*(y^a5);

end;

verbatim;

    bgp.write(M_);
    options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000000,'MaxIterations',100000,'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-8,'StepTolerance',1e-8);
    y = 1+(rand(3,1)-.5)*.5;
    g = 1+(rand(3,1)-.5)*.1;
    [y, fval, exitflag] = fsolve(@nk.bgpfun, [y;g], options);
    assert(max(abs(y-1))<1e-9);

end;