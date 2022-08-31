@#define extended_path_version = 1

var Capital, Output, Labour, Consumption,  Investment, Output1, Labour1, Consumption1, Output2, Labour2, Consumption2, Efficiency, efficiency, ExpectedTerm, LagrangeMultiplier;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, psi, delta, rho, effstar, sigma;

@#include "rbcii-calibration.inc"

model(use_dll);

  efficiency = rho*efficiency(-1) + sigma*EfficiencyInnovation;

  Efficiency = effstar*exp(efficiency);

  (((Consumption1^theta)*((1-Labour1)^(1-theta)))^(1-tau))/Consumption1 - beta*ExpectedTerm(1) + LagrangeMultiplier(1)*beta*(1-delta);

  ExpectedTerm = ((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital(-1))^(1-psi))+1-delta);

  LagrangeMultiplier = max(0, (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption - beta*ExpectedTerm(1) + LagrangeMultiplier(1)*beta*(1-delta));

  ((1-theta)/theta)*(Consumption1/(1-Labour1)) - (1-alpha)*(Output1/Labour1)^(1-psi);

  Output1 = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour1^psi))^(1/psi);

  Consumption2  = Output2;

  ((1-theta)/theta)*(Consumption2/(1-Labour2)) - (1-alpha)*(Output2/Labour2)^(1-psi);

  Output2 = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour2^psi))^(1/psi);

  Consumption = (Output1 > Consumption1)*Consumption1 + (1-(Output1 > Consumption1))*Consumption2;

  Labour = (Output1 > Consumption1)*Labour1 + (1-(Output1 > Consumption1))*Labour2;

  Output = (Output1 > Consumption1)*Output1 + (1-(Output1 > Consumption1))*Output2;

  Capital = Output-Consumption + (1-delta)*Capital(-1);

  Investment = Capital - (1-delta)*Capital(-1);

end;

@#if extended_path_version

    shocks;
    var EfficiencyInnovation = 1;
    end;

    steady(nocheck);

    options_.ep.stochastic.order = 0;
!*
    ts = extended_path([], 200, [], options_, M_, oo_);
    ts.save('rbcii-sim-data');

    options_.ep.stochastic.order = 1;
    ts1_4 = extended_path([], 200, [], options_, M_, oo_);

@#else

    shocks;
    var EfficiencyInnovation;
    periods 1;
    values -.8;
    end;

    steady;//(nocheck);

    options_.simul.maxit = 100;

    simul(periods=4000);

    n = 100;

    figure('Name','(rbcii) Investment.');
    plot(Output(1:n)-Consumption(1:n),'-b','linewidth',2)

@#endif
