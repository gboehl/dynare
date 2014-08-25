var V lL Vkp Valphaexp lzn lzd lp0 lpi MC lwreal lY lDisp lC lpiavg Int lA lDZ lG pistar Intr
lsdf Int1 lpi1;
// pricebond pricebondrn ytm ytmrn termprem ehpr slope;
varexo epsA epsZ epsG epsPistar epsInt;

parameters theta xi beta phi alpha eta KBar chi0 LMax chi IBar rhoinflavg taylrho DZBar taylpi tayly YBar rhoa rhoz rhog rhopistar piBar gssload consoldelta VAIMSS GBar;

  DZBar = 1.0025;
  eta = 2/3;
  IES = .11;
  phi = 1/IES;
  beta = .99 *DZBar^phi;
  delta = .02;

  Frisch = .28;
  CRRA = 110;

  LMax = 3;
//  chi0 = exp(wrealAIMSS -phi*CAIMSS) *(LMax-1)^chi; // normalize L in ss = 1 
  chi0 = 1/3;
  chi = 1/Frisch *(LMax-1);
  alpha = (CRRA - 1/(1/phi + 1/chi *(LMax-1))) *(1/(1-phi) + 1/(1-chi) *(LMax-1));

  theta = .2;
  xi = .78;
  
  K_Y = 10;
  YBar = K_Y^((1-eta)/eta);
  GBar = .17 *YBar;
  KBar = K_Y *YBar;
  IBar = (delta +(DZBar-1)) *KBar;
  piBar = 1.0;

  taylrho = .73;
  taylpi = .53;
  tayly = .93;

  rhoa = .96;
  rhoz = 0;
  rhog = .95;
  rhopistar = 0;
  epsafl = 1;    // these flags turn on or off each corresponding shock in the model 
  epsapermfl = 0;
  epsgfl = 1;
  epsifl = 1;
  epspistarfl = 0;
  rhoinflavg = .7;  // for price level targeting; set rhoinflavg = .99 or .999 
  gssload = 0;    // this is the GSS theta parameter in the long-run nominal risks section of the paper 

  VAIMSS = 1;
  consoldelta = 1;

model;//(use_dll);
// Value function and Euler equation
V = exp(lC)^(1-phi) /(1-phi) + chi0 *(LMax-exp(lL))^(1-chi) /(1-chi) + beta *Vkp;
Int1 = Int/10;
lpi1 = lpi/10;
exp(lC)^-phi = beta *exp(10*Int1-10*lpi1(+1)) *exp(lC(+1))^-phi *exp(lDZ(+1))^-phi
 *(V(+1) *exp(lDZ(+1))^(1-phi) /Vkp)^-alpha;
// The following two equations define the E-Z-W-K-P certainty equivalent term
// Vkp = (E_t V(+1)^(1-alpha))^(1/(1-alpha)). It takes two equations to do this because
// Perturbation AIM sets the expected value of all equations equal to zero; E_t F(variables) = 0.
// Thus; the first equation below defines Valphaexp = E_t V(+1)^(1-alpha). The second
// equation then takes the (1-alpha)th root of this expectation.
// Finally; the scaling and unscaling of Valphaexp by the constant VAIMSS and DZBar improves the
// numerical behavior of the model; without it; the steady-state value of Valphaexp can be minuscule
// (e.g.; 10^-50); which requires Mathematica to use astronomical levels of precision to solve. *)
[mcp = 'Valphaexp > 1e-5']
Valphaexp = (V(+1) *exp(lDZ(+1))^(1-phi) /VAIMSS /DZBar^(1-phi))^(1-alpha);
Vkp = VAIMSS *DZBar^(1-phi) *Valphaexp^(1/(1-alpha));

exp(lzn) = (1+theta) *MC *exp(lY) + xi *beta *exp(lC(+1)-lC)^-phi *exp(lDZ(+1))^-phi
 *(V(+1) *exp(lDZ(+1))^(1-phi) /Vkp)^-alpha *exp(lpi(+1))^((1+theta)/theta/eta) *exp(lzn(+1));
exp(lzd) = exp(lY) + xi *beta *exp(lC(+1)-lC)^-phi*exp(lDZ(+1))^-phi
 *(V(+1) *exp(lDZ(+1))^(1-phi) /Vkp)^-alpha *exp(lpi(+1))^(1/theta) *exp(lzd(+1));
exp(lp0)^(1+(1+theta)/theta *(1-eta)/eta) = exp(lzn-lzd);
exp(lpi)^(-1/theta) = (1-xi) *exp(lp0+lpi)^(-1/theta) + xi;
// Marginal cost and real wage 
MC = exp(lwreal) /eta *exp(lY)^((1-eta)/eta) /exp(lA)^(1/eta) /KBar^((1-eta)/eta);
[mcp = 'lL < 1.0986']
chi0 *(LMax-exp(lL))^-chi /exp(lC)^-phi = exp(lwreal);
// Output equations 
exp(lY) = exp(lA) *KBar^(1-eta) *exp(lL)^eta /exp(lDisp);
exp(lDisp)^(1/eta) = (1-xi) *exp(lp0)^(-(1+theta)/theta/eta)
 + xi *exp(lpi)^((1+theta)/theta/eta) *exp(lDisp(-1))^(1/eta);
exp(lC) = exp(lY)-exp(lG)-IBar; // aggregate resource constraint; no adj costs 
// Monetary Policy Rule 
lpiavg = rhoinflavg *lpiavg(-1) + (1-rhoinflavg) *lpi;
4*Int = (1-taylrho) * ( 4*log(1/beta *DZBar^phi) + 4*lpiavg
 + taylpi * (4*lpiavg-pistar) + tayly * (exp(lY)-YBar)/YBar )
 + taylrho * 4*Int(-1) + epsInt; // multiply Int; infl by 4 to put at annual rate 
// Exogenous Shocks 
lA = rhoa * lA(-1) + epsA;
lDZ = (1-rhoz)*log(DZBar) + rhoz * lDZ(-1) + epsZ;
lG = (1-rhog)*log(GBar) + rhog * lG(-1) + epsG;
pistar = (1-rhopistar) *log(piBar) + rhopistar *pistar(-1) + gssload *(4*lpiavg-pistar)
 + epsPistar;
// Term premium and other auxiliary finance equations 
Intr = Int(-1) - lpi; // ex post real short rate 
exp(lsdf) = beta *exp(lC(+1)-lC)^-phi *exp(lDZ(+1))^-phi
 *(V(+1) *exp(lDZ(+1))^(1-phi) /Vkp)^-alpha /exp(lpi(+1));
//pricebond = 1 + consoldelta *beta *exp(lC(+1)-lC)^-phi *exp(lDZ(+1))^-phi
// *(V(+1) *exp(lDZ(+1))^(1-phi) /Vkp)^-alpha /pi(+1) *pricebond(+1);
//pricebondrn = 1 + consoldelta *pricebondrn(+1) /exp(Int);
//ytm = log(consoldelta*pricebond/(pricebond-1)) *400; // yield in annualized pct 
//ytmrn = log(consoldelta*pricebondrn/(pricebondrn-1)) *400;
//termprem = 100 * (ytm–ytmrn); // term prem in annualized basis points 
//ehpr = ( (consoldelta *pricebond + exp(Int(-1))) /pricebond(-1)–exp(Int(-1)) *400;
//slope = ytm–Int*400;
end;

steady_state_model;
lA = 0;
lDZ = log(DZBar);
lG = log(GBar);
pistar = log(piBar);
lpi = log(piBar);
lpiavg = lpi;
lDisp = 0;
Int = log(DZBar^phi/beta)+lpi;
Int1 = Int/10;
lpi1 = lpi/10;
lL = 0;
lY = log(YBar);
lC = log(exp(lY)-exp(lG)-IBar);
lp0 = log(((exp(lpi)^(-1/theta)-xi)/((1-xi)*exp(lpi)^(-1/theta)))^(-theta));
lzd = log(exp(lY)/(1-xi *beta *exp(lDZ)^-phi
  *exp(lpi)^(1/theta)));
lzn = log(exp(lp0)^(1+(1+theta)/theta *(1-eta)/eta)*exp(lzd));
MC = (exp(lzn)-xi *beta  *DZBar^-phi*piBar^((1+theta)/theta/eta)*exp(lzn))/((1+theta)*YBar);
lwreal = log(MC*eta /(exp(lY)^((1-eta)/eta) / KBar^((1-eta)/eta)));

chi0 = exp(lwreal)*exp(lC)^(-phi)/(LMax-exp(lL))^(-chi);
V = (exp(lC)^(1-phi) /(1-phi) + chi0 *(LMax-exp(lL))^(1-chi) /(1-chi))/(1-beta*DZBar^(1-phi));
k = 1.05;
VAIMSS = k*V;
Valphaexp = k^(alpha-1);
Vkp = DZBar^(1-phi)*V;
Intr = Int - lpi;
lsdf = log(beta *exp(lDZ)^-phi/exp(lpi));
end;

steady;

shocks;
var epsA; stderr 0.005;
var epsZ; stderr 0.001;
var epsG; stderr 0.004;
var epsPistar; stderr 0.0005;
var epsInt; stderr 0.003;
end;

//stoch_simul(order=3,periods=50000,pruning);

options_.ep.verbosity=0;
options_.ep.stochastic.algo=1;
options_.ep.solve_algo = 10;
options_.ep.maxit = 100;
options_.ep.IntegrationAlgorithm='UT_2p+1';
options_.ep.ut.k = 1;
options_.solve_tolf = 1e-12;

extended_path(order=1,periods=3);