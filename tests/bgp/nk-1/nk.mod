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
    if isoctave
        options = optimset('Display', 'iter', 'MaxFunEvals', 1000000,'MaxIter',100000,'Jacobian','on','TolFun',1e-8,'TolX',1e-8);
    else
        options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000000,'MaxIterations',100000,'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-8,'StepTolerance',1e-8);
    end
    if isoctave && octave_ver_less_than('6')
        % Octave < 6 can't take a function handle of a function in a package
        % See https://savannah.gnu.org/bugs/index.php?46659
        fun = str2func('nk.bgpfun');
    else
        fun = @nk.bgpfun;
    end
    y = 1+(rand(3,1)-.5)*.5;
    g = 1+(rand(3,1)-.5)*.1;
    [y, fval, exitflag] = fsolve(fun, [y;g], options);
    assert(max(abs(y-1))<1e-9);

end;
