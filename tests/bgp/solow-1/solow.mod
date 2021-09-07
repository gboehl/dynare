// This mod file test the bgp.write function by characterizing the Balanced Growth Path of the Solow model. This is done 10000 times in a Monte Carlo loop
// randomizing the initial guess of the nonlinear solver.

var Efficiency
    EfficiencyGrowth
    Population
    PopulationGrowth
    Output
    PhysicalCapitalStock ;

varexo e_x
       e_n ;

parameters alpha
	       delta
	       s
           rho_x
           rho_n
           EfficiencyGrowth_ss
           PopulationGrowth_ss ;

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.01;
PopulationGrowth_ss = 1.01;

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = PhysicalCapitalStock(-1)^alpha*(Efficiency*Population)^(1-alpha);
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

if ~isoctave
% The levenberg-marquardt algorithm is not available in octave's
% implementation of fsolve, so we skip this check for Octave
verbatim;

    bgp.write(M_);
    MC = 10000;
    KY = NaN(MC,1);
    GY = NaN(MC,1);
    GK = NaN(MC,1);
    EG = NaN(MC,1);
    if isoctave
        options = optimset('Display', 'off', 'MaxFunEvals', 1000000,'MaxIter',100000,'Jacobian','on','TolFun',1e-8,'TolX',1e-8);
    elseif matlab_ver_less_than('9.0')
        % See https://fr.mathworks.com/help/optim/ug/current-and-legacy-option-name-tables.html
        options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt','MaxFunEvals',1000000,'MaxIter',100000,'Jacobian','on','TolFun',1e-8,'TolX',1e-8);
    else
        options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000000,'MaxIterations',100000,'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-8,'StepTolerance',1e-8);
    end
    if isoctave && octave_ver_less_than('6')
        % Octave < 6 can't take a function handle of a function in a package
        % See https://savannah.gnu.org/bugs/index.php?46659
        fun = str2func('solow.bgpfun');
    else
        fun = @solow.bgpfun;
    end
    reverseStr = '';
    for i=1:MC
        y = 1+(rand(6,1)-.5)*.2;
        g = ones(6,1);
        [y, fval, exitflag] = fsolve(fun, [y;g], options);
        if exitflag>0
            KY(i) = y(6)/y(5);
            GY(i) = y(11);
            GK(i) = y(12);
            EG(i) = y(2);
        end        
        % Display the progress
        percentDone = 100 * i / MC;
        msg = sprintf('Percent done: %3.1f', percentDone);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
    % Compute the physical capital stock over output ratio along the BGP as
    % a function of the deep parameters...
    theoretical_long_run_ratio = s*EfficiencyGrowth_ss*PopulationGrowth_ss/(EfficiencyGrowth_ss*PopulationGrowth_ss-1+delta);
    % ... and compare to the average ratio in the Monte Carlo
    mu = mean(KY(~isnan(KY)));
    s2 = var(KY(~isnan(KY)));
    nn = sum(isnan(KY))/MC;
    disp(sprintf('Theoretical K/Y BGP ratio = %s.', num2str(theoretical_long_run_ratio)));
    disp(sprintf('mean(K/Y) = %s.', num2str(mu)));
    disp(sprintf('var(K/Y) = %s.', num2str(s2)));
    disp(sprintf('Number of failures: %u (over %u problems).', sum(isnan(KY)), MC));
    assert(abs(mu-theoretical_long_run_ratio)<1e-8);
    assert(s2<1e-16);
    assert(abs(mean(GY(~isnan(GY)))-EfficiencyGrowth_ss*PopulationGrowth_ss)<1e-8)
    assert(var(GY(~isnan(GY)))<1e-16);
    assert(abs(mean(GK(~isnan(GK)))-EfficiencyGrowth_ss*PopulationGrowth_ss)<1e-8)
    assert(var(GK(~isnan(GK)))<1e-16);
    assert(abs(mean(EG(~isnan(EG)))-EfficiencyGrowth_ss)<1e-8)
    assert(var(EG(~isnan(EG)))<1e-16);

end;

end;
