/*
 * Model codes kindly provided by Mathias Trabandt, modified to work with
 * Parameters here are set such that the AOB Model is active.
 *
 * Please note that the following copyright notice only applies to this Dynare
 * implementation of the model.
 */

/*
 * Copyright © 2023 Dynare Team
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <https://www.gnu.org/licenses/>.
 */


/*
 * Wanted: monetary policy shock is not in period t information set of agents 
 *         (just as in VAR identified with Choleski decomp.) 
 *         All other shocks are in period t information set. 
 * Solution: create two parallel economies which are structurally identical.
 *           One which is subject to the monetary shock only.
 *           This shock is not in period t info set.
 *           The other economy has the other shocks with the standard info set.
 * Naming convention: variables ending with 'F' are in the economy with full info set
 *                    variables ending with 'R' are in the restricted info set economy
 *                    variables ending with 'AGG' are aggregated variables from both economies
 * Note that in the 'R' economy, all exogenous variables are set to their steady states, except the monetary shock
 * Likewise, in the 'F' economy, the monetary shock is set to zero while all other shocks are active
 * Both steady states (i.e. 'R' and 'F' economy) are identical
 * All variables are in logs
 */

// Endogenous variables (R economy)
var
wpR wpF wR wF psiR cR RR piR iR pkprimeR FR KR ukR
kbarR lR xR varthetR phaloR yR JR UR GDPR fR vR
unempR VR AAR vTotR QR piwR whaloR FwR KwR 
varthetpR
;  

// Endogenous variables (F economy)
var
psiF cF RF piF iF pkprimeF FF KF ukF
kbarF lF xF varthetF phaloF yF JF UF GDPF fF 
vF unempF VF AAF vTotF QF piwF whaloF FwF KwF
muzF mupsiF muF nGF nPHIF nGAMF nDF nKAPF 
varthetpF
;

// Endogenous variables (AGG economy)
var
GDPAGG piAGG RAGG ukAGG lAGG wAGG cAGG 
iAGG unempAGG vTotAGG fAGG pinvestAGG
;

// Shocks
varexo
epsR_eps muz_eps mupsi_eps
;

parameters
ydiffdataPercent,idiffdataPercent,alfa,rho,u,Q,sigm,
betta,lambda,deltak,tau,nuf,etag,kappa,b,xi,D,delta,gamma,Spp 
sigmab sigmaa phi rhoR rpi ry sig_epsR sigmam kappaf DSHARE deltapercent
recSHAREpercent iota doNash eta doEHL xiw lambdaw AEHL  
f rhomupsi rhomuz sig_mupsi sig_muz thetaG thetaPHI thetaKAP thetaD 
dolagZBAR profy varkappaw pibreve thetaw varkappaf M sigmaL 
tau_SS epsilon_SS upsilon_SS zetac_SS g_SS pibar thetaGAM kappaw
unemp_SS vTot_SS alp1 alp2 alp3 alp4 bet1 bet2 bet3 s searchSHAREpercent
;
   
model;

// auxiliary equations for R economy such that current realization fo monetary policy shock is not in information set
//information set in line with Choleski decomposition in VAR
//set exogenous variables to steady states (except monetary shock, of course)
#epsilonR=epsilon_SS;
#upsilonR=upsilon_SS;
#upsilonR_tp1=upsilon_SS;
#zetacR=zetac_SS;
#zetacR_tp1=zetac_SS;
#gR=g_SS;
#tauR=tau_SS;
#mupsiR=STEADY_STATE(mupsiF);
#mupsiR_tp1=STEADY_STATE(mupsiF);
#muzR=STEADY_STATE(muzF);
#muR=alfa/(1-alfa)*STEADY_STATE(mupsiF)+STEADY_STATE(muzF); 
#muR_tp1=STEADY_STATE(muF);
#nGR=STEADY_STATE(nGF);
#nPHIR=STEADY_STATE(nPHIF);
#nGAMR=STEADY_STATE(nGAMF);
#nDR=STEADY_STATE(nDF);
#nKAPR=STEADY_STATE(nKAPF);

//conditional expectations based on t-1 info;
#piR_tp1=EXPECTATION(-1)(piR(+1));
#FR_tp1=EXPECTATION(-1)(FR(+1));
#KR_tp1=EXPECTATION(-1)(KR(+1));
#cR_tp1=EXPECTATION(-1)(cR(+1));
#psiR_tp1=EXPECTATION(-1)(psiR(+1));
#pkprimeR_tp1=EXPECTATION(-1)(pkprimeR(+1));
#iR_tp1=EXPECTATION(-1)(iR(+1));
#UR_tp1=EXPECTATION(-1)(UR(+1));
#fR_tp1=EXPECTATION(-1)(fR(+1));
#wR_tp1=EXPECTATION(-1)(wR(+1));
#VR_tp1=EXPECTATION(-1)(VR(+1)); 
#ukR_tp1=EXPECTATION(-1)(ukR(+1));
#piwR_tp1=EXPECTATION(-1)(piwR(+1));
#KwR_tp1=EXPECTATION(-1)(KwR(+1));
#FwR_tp1=EXPECTATION(-1)(FwR(+1));
#wpR_tp1=EXPECTATION(-1)(wpR(+1));
#varthetpR_tp1=EXPECTATION(-1)(varthetpR(+1));
#AAR_tp1=EXPECTATION(-1)(AAR(+1));
#RRinfo=EXPECTATION(-1)(RR);

// abbreviations
#aofukprimeR=sigmab*sigmaa*(exp(ukR))+sigmab*(1-sigmaa);
#aofukprimeR_tp1=sigmab*sigmaa*(exp(ukR_tp1))+sigmab*(1-sigmaa);
#aofukR=0.5*sigmab*sigmaa*(exp(ukR))^2+sigmab*(1-sigmaa)*exp(ukR)+sigmab*((sigmaa/2)-1);
#aofukR_tp1=0.5*sigmab*sigmaa*(exp(ukR_tp1))^2+sigmab*(1-sigmaa)*exp(ukR_tp1)+sigmab*((sigmaa/2)-1);
#RkR_tp1=log(((exp(ukR_tp1)*aofukprimeR_tp1-aofukR_tp1)+(1-deltak)*exp(pkprimeR_tp1))/(exp(pkprimeR+mupsiR_tp1-piR_tp1)));
#StildeR=(0.5*(exp(sqrt(Spp)*(exp(iR)/exp(iR(-1))*exp(muR)*exp(mupsiR)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))+exp(-sqrt(Spp)*(exp(iR)/exp(iR(-1))*exp(muR)*exp(mupsiR)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-2));
#StildeprimeR=(0.5*sqrt(Spp)*(exp(sqrt(Spp)*(exp(iR)/exp(iR(-1))*exp(muR)*exp(mupsiR)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-exp(-sqrt(Spp)*(exp(iR)/exp(iR(-1))*exp(muR)*exp(mupsiR)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))));
#StildeprimeR_tp1=(0.5*sqrt(Spp)*(exp(sqrt(Spp)*(exp(iR_tp1)/exp(iR)*exp(muR_tp1)*exp(mupsiR_tp1)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-exp(-sqrt(Spp)*(exp(iR_tp1)/exp(iR)*exp(muR_tp1)*exp(mupsiR_tp1)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))));
#mcR=tauR+alfa*mupsiR+( (1-doEHL)*varthetR+doEHL*wR) +log(nuf*(exp(RRinfo))+1-nuf)-epsilonR-log(1-alfa)-alfa*(kbarR(-1)+ukR-muR-lR-(doEHL*(lambdaw/(lambdaw-1))*whaloR) ); 
#pitildewpiR=-piwR+kappaw*piR(-1)+(1-kappaw-varkappaw)*log(pibar)+varkappaw*log(pibreve)+thetaw*STEADY_STATE(muF);
#pitildewpiR_tp1=-piwR_tp1+kappaw*piR+(1-kappaw-varkappaw)*log(pibar)+varkappaw*log(pibreve)+thetaw*STEADY_STATE(muF);
#pitildepiR=-piR+kappaf*piR(-1)+(1-kappaf-varkappaf)*log(pibar)+varkappaf*log(pibreve);
#pitildepiR_tp1=-piR_tp1+kappaf*piR+(1-kappaf-varkappaf)*log(pibar)+varkappaf*log(pibreve);

//R1 - Consumption FOC
exp(psiR)=exp(zetacR)/(exp(cR)-b*exp(cR(-1))/exp(muR))-betta*b*(exp(zetacR_tp1)/
(exp(cR_tp1)*exp(muR_tp1)-b*exp(cR)));

//R2 - Bond FOC (Fisher equation) 
exp(psiR)=betta*exp(psiR_tp1)/exp(muR_tp1)*exp(RRinfo)/exp(piR_tp1);

//R3 - Investment FOC
1=exp(pkprimeR)*exp(upsilonR)*(1-StildeR-StildeprimeR*exp(muR)*exp(mupsiR)*exp(iR)/exp(iR(-1)))+betta*exp(psiR_tp1)/exp(psiR)*exp(pkprimeR_tp1)*exp(upsilonR_tp1)*StildeprimeR_tp1*(exp(iR_tp1)/exp(iR))^2*exp(muR_tp1)*exp(mupsiR_tp1);

//R4 - Capital FOC
1=betta*exp(psiR_tp1-psiR+RkR_tp1-muR_tp1-piR_tp1);

//R5 - Law of motion for physical capital
exp(kbarR)=(1-deltak)/(exp(muR)*exp(mupsiR))*exp(kbarR(-1))+exp(upsilonR)*(1-StildeR)*exp(iR);

//R6 - Cost minimization (opt. factor inputs)
0=aofukprimeR*exp(ukR)*exp(kbarR(-1))/(exp(muR)*exp(mupsiR))-alfa/(1-alfa)*(nuf*exp(RRinfo)+1-nuf)*exp(lR)*((1-doEHL)*(exp(varthetR))+doEHL*( exp(wR)*(exp(whaloR))^(lambdaw/(lambdaw-1))));

//R7 - Production
exp(yR)=exp(phaloR)^(lambda/(lambda-1))*(exp(epsilonR)*((doEHL*exp(whaloR)^(lambdaw/(lambdaw-1))+(1-doEHL)*1)*exp(lR))^(1-alfa)*
             (exp(kbarR(-1)+ukR)/(exp(muR)*exp(mupsiR)))^alfa-phi*exp(nPHIR)); 

//R8 - Resource Constraint
exp(yR)=exp(gR)*exp(nGR)+exp(cR)+exp(iR)+aofukR*exp(kbarR(-1))/(exp(mupsiR)*exp(muR)) 
  +(1-doEHL)*( s*exp(nKAPR)*exp(xR)/exp(QR)*exp(lR(-1))  +  kappa*exp(nKAPR)*exp(xR)*exp(lR(-1))   );

//R9 Monetary Policy Rule
RR=rhoR*RR(-1)+(1-rhoR)*(STEADY_STATE(RR)+rpi*(piR-STEADY_STATE(piF))+ry*(GDPR-STEADY_STATE(GDPF)))-sig_epsR*epsR_eps/400;

//R10 Pricing 1
exp(FR)-exp(psiR+yR)-betta*xi*exp( pitildepiR_tp1  /(1-lambda)+FR_tp1);

//R11 Pricing 2 
exp(KR)-lambda*exp(psiR+yR+mcR)-betta*xi*exp( pitildepiR_tp1 *lambda/(1-lambda)+KR_tp1);

//R12 Pricing 3
KR-FR=(1-lambda)*log((1-xi*exp( pitildepiR /(1-lambda)))/(1-xi));

//R13  Price dispersion
exp(phaloR*lambda/(1-lambda))-(1-xi)^(1-lambda)*(1-xi*exp( pitildepiR /(1-lambda)))^lambda
   -xi*(exp( pitildepiR +phaloR(-1)))^(lambda/(1-lambda));

//R14 Present value of wages
doEHL*(exp(wpR)-exp(STEADY_STATE(wpF)))=(1-doEHL)*(-exp(wpR)+exp(wR)+rho*betta*exp(psiR_tp1)/exp(psiR)*exp(wpR_tp1));

//R15 Present value of marginal revenue
doEHL*(exp(varthetpR)-exp(STEADY_STATE(varthetpF)))=(1-doEHL)*(-exp(varthetpR)+exp(varthetR)+rho*betta*exp(psiR_tp1)/exp(psiR)*exp(varthetpR_tp1));

//R16 Hiring FOC  - zero profit/free entry condition
doEHL*(exp(xR)-exp(STEADY_STATE(xF)))=(1-doEHL)*(  -s*exp(nKAPR)/exp(QR) -kappa*exp(nKAPR) + exp(JR)   );

//R17 Value of firm
doEHL*(exp(JR)-exp(STEADY_STATE(JF)))=(1-doEHL)*(-exp(JR)+exp(varthetpR)-exp(wpR));

//R18 Value of work
doEHL*(exp(VR)-exp(STEADY_STATE(VF)))=(1-doEHL)*(-exp(VR)+exp(wpR)+exp(AAR));

//R19 Present value of worker payoff after separation 
doEHL*(exp(AAR)-exp(STEADY_STATE(AAF)))=(1-doEHL)*(-exp(AAR)+(1-rho)*betta*exp(psiR_tp1)/exp(psiR)*(exp(fR_tp1)*exp(VR_tp1)
+(1-exp(fR_tp1))*exp(UR_tp1))+rho*betta*exp(psiR_tp1)/exp(psiR)*exp(AAR_tp1));

//R20 Unemployment value
doEHL*(exp(UR)-exp(STEADY_STATE(UF)))=(1-doEHL)*(-exp(UR)+D*exp(nDR)+betta*exp(psiR_tp1)/exp(psiR)*(exp(fR_tp1)*exp(VR_tp1)+(1-exp(fR_tp1))*exp(UR_tp1)));

//R21 Sharing rule
doEHL*(exp(varthetR)-STEADY_STATE(exp(varthetF)))=(1-doEHL)*(doNash*(exp(VR)-exp(UR)-eta*(exp(JR)+exp(VR)-exp(UR)))+(1-doNash)*(
  
-exp(JR)+bet1*(exp(VR)-exp(UR))-bet2*exp(nGAMR)*gamma+bet3*(exp(varthetR)-exp(nDR)*D)

));

//R22 GDP
exp(GDPR)=exp(gR)*exp(nGR)+exp(cR)+exp(iR);

//R23 Unempl. rate
doEHL*(exp(unempR)-exp(STEADY_STATE(unempF)))=(1-doEHL)*(-exp(unempR)+1-exp(lR));

//R24 Job finding rate
doEHL*(exp(fR)-exp(STEADY_STATE(fF)))=(1-doEHL)*(-exp(fR)+exp(xR)*exp(lR(-1))/(1-rho*exp(lR(-1)))); 

//R25 Matching function
doEHL*(exp(vTotR)-exp(STEADY_STATE(vTotF)))=(1-doEHL)*(-exp(xR)*exp(lR(-1))+sigmam*(1-rho*exp(lR(-1)))^sigm*(exp(vTotR))^(1-sigm));

//R26 Total vacancies
doEHL*(exp(vR)-exp(STEADY_STATE(vF)))=(1-doEHL)*(-exp(vTotR)+exp(vR)*exp(lR(-1)));

//R27 Vacancy filling rate
doEHL*(exp(QR)-exp(STEADY_STATE(QF)))=(1-doEHL)*(-exp(QR)+exp(xR)/exp(vR));

//R28 Wage inflation (EHL)
exp(piwR)=exp(piR)*exp(muR)*exp(wR)/exp(wR(-1));

//R29  Wage dispersion
doEHL*(exp(whaloR*lambdaw/(1-lambdaw))-(1-xiw)^(1-lambdaw)*(1-xiw*exp(pitildewpiR/(1-lambdaw)))^lambdaw
   -xiw*(exp(pitildewpiR+whaloR(-1)))^(lambdaw/(1-lambdaw)))=(1-doEHL)*(exp(whaloR)-exp(STEADY_STATE(whaloF)));

//R30 wage setting 1 (EHL)
doEHL*(-exp(FwR) + exp(psiR)/lambdaw*exp(whaloR)^( lambdaw/(lambdaw-1))*exp(lR) + betta*xiw*(exp(wR_tp1)/exp(wR))*(exp(pitildewpiR_tp1))^( 1/(1-lambdaw) )*exp(FwR_tp1))
=(1-doEHL)*(exp(FwR)-exp(STEADY_STATE(FwF)));

//R31 Law of motion of employment or EHL wage setting 2
doEHL*(-exp(KwR) + ( exp(whaloR)^( lambdaw/(lambdaw-1))*exp(lR))^(1+sigmaL) + betta*xiw*(exp(pitildewpiR_tp1))^( lambdaw*(1+sigmaL)/(1-lambdaw) )*exp(KwR_tp1))    
=(1-doEHL)*(-exp(lR)+(rho+exp(xR))*exp(lR(-1)));

//R32 wage setting 3 (EHL)
doEHL*(1-(1-xiw)*(AEHL*exp(KwR)/(exp(FwR)*exp(wR)))^( 1/( 1-lambdaw*(1+sigmaL) )) - xiw*(exp(pitildewpiR))^(1/(1-lambdaw)))=(1-doEHL)*(exp(KwR)-exp(STEADY_STATE(KwF)));

////////////////////////////////////////////////////////////////////// 
//F economy, other shocks, current realization in info set /////////// 
////////////////////////////////////////////////////////////////////// 
//set not needed exogenous variables to steady states 
#tauF=tau_SS;
#upsilonF=upsilon_SS;
#upsilonF_tp1=upsilon_SS;
#epsilonF=epsilon_SS;
#zetacF=zetac_SS;
#zetacF_tp1=zetac_SS;
#gF=g_SS;

//abbreviations
#aofukprimeF=sigmab*sigmaa*(exp(ukF))+sigmab*(1-sigmaa);
#aofukprimeF_tp1=sigmab*sigmaa*(exp(ukF(+1)))+sigmab*(1-sigmaa);
#aofukF=0.5*sigmab*sigmaa*(exp(ukF))^2+sigmab*(1-sigmaa)*exp(ukF)+sigmab*((sigmaa/2)-1);
#aofukF_tp1=0.5*sigmab*sigmaa*(exp(ukF(+1)))^2+sigmab*(1-sigmaa)*exp(ukF(+1))+sigmab*((sigmaa/2)-1);
#RkF_tp1=log(((exp(ukF(+1))*aofukprimeF_tp1-aofukF_tp1)+(1-deltak)*exp(pkprimeF(+1)))/(exp(pkprimeF+mupsiF(+1)-piF(+1))));
#StildeF=(0.5*(exp(sqrt(Spp)*(exp(iF)/exp(iF(-1))*exp(muF)*exp(mupsiF)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))+exp(-sqrt(Spp)*(exp(iF)/exp(iF(-1))*exp(muF)*exp(mupsiF)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-2));
#StildeprimeF=(0.5*sqrt(Spp)*(exp(sqrt(Spp)*(exp(iF)/exp(iF(-1))*exp(muF)*exp(mupsiF)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-exp(-sqrt(Spp)*(exp(iF)/exp(iF(-1))*exp(muF)*exp(mupsiF)-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))));
#StildeprimeF_tp1=(0.5*sqrt(Spp)*(exp(sqrt(Spp)*(exp(iF(+1))/exp(iF)*exp(muF(+1))*exp(mupsiF(+1))-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))-exp(-sqrt(Spp)*(exp(iF(+1))/exp(iF)*exp(muF(+1))*exp(mupsiF(+1))-exp(STEADY_STATE(muF))*exp(STEADY_STATE(mupsiF))))));
#mcF=tauF+alfa*mupsiF+( (1-doEHL)*varthetF+doEHL*wF) +log(nuf*(exp(RF))+1-nuf)-epsilonF-log(1-alfa)-alfa*(kbarF(-1)+ukF-muF-lF-(doEHL*(lambdaw/(lambdaw-1))*whaloF) ); 
#pitildewpiF=-piwF+kappaw*piF(-1)+(1-kappaw-varkappaw)*log(pibar)+varkappaw*log(pibreve)+thetaw*STEADY_STATE(muF);
#pitildewpiF_tp1=-piwF(+1)+kappaw*piF+(1-kappaw-varkappaw)*log(pibar)+varkappaw*log(pibreve)+thetaw*STEADY_STATE(muF);
#pitildepiF=-piF+kappaf*piF(-1)+(1-kappaf-varkappaf)*log(pibar)+varkappaf*log(pibreve);
#pitildepiF_tp1=-piF(+1)+kappaf*piF+(1-kappaf-varkappaf)*log(pibar)+varkappaf*log(pibreve);

//F1 - Consumption FOC
exp(psiF)=exp(zetacF)/(exp(cF)-b*exp(cF(-1))/exp(muF))-betta*b*(exp(zetacF_tp1)/
(exp(cF(+1))*exp(muF(+1))-b*exp(cF)));

//F2 - Bond FOC (Fisher equation) 
exp(psiF)=betta*exp(psiF(+1))/exp(muF(+1))*exp(RF)/exp(piF(+1));

//F3 - Investment FOC
1=exp(pkprimeF)*exp(upsilonF)*(1-StildeF-StildeprimeF*exp(muF)*exp(mupsiF)*exp(iF)/exp(iF(-1)))+betta*exp(psiF(+1))/exp(psiF)*exp(pkprimeF(+1))*exp(upsilonF_tp1)*StildeprimeF_tp1*(exp(iF(+1))/exp(iF))^2*exp(muF(+1))*exp(mupsiF(+1));

//F4 - Capital FOC
1=betta*exp(psiF(+1)-psiF+RkF_tp1-muF(+1)-piF(+1));

//F5 - Law of motion for physical capital
exp(kbarF)=(1-deltak)/(exp(muF)*exp(mupsiF))*exp(kbarF(-1))+exp(upsilonF)*(1-StildeF)*exp(iF);

//F6 - Cost minimization (opt. factor inputs)
0=aofukprimeF*exp(ukF)*exp(kbarF(-1))/(exp(muF)*exp(mupsiF))-alfa/(1-alfa)*(nuf*exp(RF)+1-nuf)*exp(lF)*((1-doEHL)*(exp(varthetF))+doEHL*( exp(wF)*(exp(whaloF))^(lambdaw/(lambdaw-1))));

//F7 - Production
exp(yF)=exp(phaloF)^(lambda/(lambda-1))*(exp(epsilonF)*(     (doEHL*exp(whaloF)^(lambdaw/(lambdaw-1))+(1-doEHL)*1)*exp(lF)     )^(1-alfa)*
             (exp(kbarF(-1)+ukF)/(exp(muF)*exp(mupsiF)))^alfa-phi*exp(nPHIF)); 

//F8 - Resource Constraint
exp(yF)=exp(gF)*exp(nGF)+exp(cF)+exp(iF)+aofukF*exp(kbarF(-1))/(exp(mupsiF)*exp(muF)) 
  +(1-doEHL)*( s*exp(nKAPF)*exp(xF)/exp(QF)*exp(lF(-1))  +  kappa*exp(nKAPF)*exp(xF)*exp(lF(-1))   );


//F9 Monetary Policy Rule
RF=rhoR*RF(-1)+(1-rhoR)*(STEADY_STATE(RF)+rpi*(piF-STEADY_STATE(piF))+ry*(GDPF-STEADY_STATE(GDPF)));

//F10 Pricing 1
exp(FF)-exp(psiF+yF)-betta*xi*exp( pitildepiF_tp1 /(1-lambda)+FF(+1));

//F11 Pricing 2 
exp(KF)-lambda*exp(psiF+yF+mcF)-betta*xi*exp( pitildepiF_tp1*lambda/(1-lambda)+KF(+1));

//F12 Pricing 3
KF-FF=(1-lambda)*log((1-xi*exp( pitildepiF /(1-lambda)))/(1-xi));

//F13  Price dispersion
exp(phaloF*lambda/(1-lambda))-(1-xi)^(1-lambda)*(1-xi*exp( pitildepiF /(1-lambda)))^lambda
   -xi*(exp( pitildepiF +phaloF(-1)))^(lambda/(1-lambda));

//F14 Present value of wages
doEHL*(exp(wpF)-STEADY_STATE(exp(wpF)))=(1-doEHL)*(-exp(wpF)+exp(wF)+rho*betta*exp(psiF(+1))/exp(psiF)*exp(wpF(+1)));

//F15 Present value of marginal revenue
doEHL*(exp(varthetpF)-STEADY_STATE(exp(varthetpF)))=(1-doEHL)*(-exp(varthetpF)+exp(varthetF)+rho*betta*exp(psiF(+1))/exp(psiF)*exp(varthetpF(+1)));

//F16 Hiring FOC  - zero profit/free entry condition
doEHL*(exp(xF)-exp(STEADY_STATE(xF)))=(1-doEHL)*(  -s*exp(nKAPF)/exp(QF) - kappa*exp(nKAPF) + exp(JF)   );

//F17 Firm surplus
doEHL*(exp(JF)-STEADY_STATE(exp(JF)))=(1-doEHL)*(-exp(JF)+exp(varthetpF)-exp(wpF));

//F18 Value of work
doEHL*(exp(VF)-exp(STEADY_STATE(VF)))=(1-doEHL)*(-exp(VF)+exp(wpF)+exp(AAF));

//F19 Present value of worker payoff after separation 
doEHL*(exp(AAF)-exp(STEADY_STATE(AAF)))=(1-doEHL)*(-exp(AAF)+(1-rho)*betta*exp(psiF(+1))/exp(psiF)*(exp(fF(+1))*exp(VF(+1))+(1-exp(fF(+1)))*exp(UF(+1)))+rho*betta*exp(psiF(+1))/exp(psiF)*exp(AAF(+1)));
 
//F20 Unempoyment value
doEHL*(exp(UF)-exp(STEADY_STATE(UF)))=(1-doEHL)*(-exp(UF)+D*exp(nDF)+betta*exp(psiF(+1))/exp(psiF)*(exp(fF(+1))*exp(VF(+1))+(1-exp(fF(+1)))*exp(UF(+1))));

//R21 Sharing rule
doEHL*(exp(varthetF)-STEADY_STATE(exp(varthetF)))=(1-doEHL)*(doNash*(exp(VF)-exp(UF)-eta*(exp(JF)+exp(VF)-exp(UF)))+(1-doNash)*(

-exp(JF)+bet1*(exp(VF)-exp(UF))-bet2*exp(nGAMF)*gamma+bet3*(exp(varthetF)-exp(nDF)*D)

));

//F22 GDP
exp(GDPF)=exp(gF)*exp(nGF)+exp(cF)+exp(iF);

//F23 Unempl. rate
doEHL*(exp(unempF)-exp(STEADY_STATE(unempF)))=(1-doEHL)*(-exp(unempF)+1-exp(lF));

//F24 Finding rate
doEHL*(exp(fF)-exp(STEADY_STATE(fF)))=(1-doEHL)*(-exp(fF)+exp(xF)*exp(lF(-1))/(1-rho*exp(lF(-1)))); 

//F25 Matching function
doEHL*(exp(vTotF)-exp(STEADY_STATE(vTotF)))=(1-doEHL)*(-exp(xF)*exp(lF(-1))+sigmam*(1-rho*exp(lF(-1)))^sigm*(exp(vTotF))^(1-sigm));

//F26 Total vacancies
doEHL*(exp(vF)-exp(STEADY_STATE(vF)))=(1-doEHL)*(-exp(vTotF)+exp(vF)*exp(lF(-1)));

//F27 Vacancy filling rate
doEHL*(exp(QF)-exp(STEADY_STATE(QF)))=(1-doEHL)*(-exp(QF)+exp(xF)/exp(vF));

//F28 Wage inflation
exp(piwF)=exp(piF)*exp(muF)*exp(wF)/exp(wF(-1));

//F29  Wage dispersion
doEHL*(exp(whaloF*lambdaw/(1-lambdaw))-(1-xiw)^(1-lambdaw)*(1-xiw*exp(pitildewpiF/(1-lambdaw)))^lambdaw
   -xiw*(exp(pitildewpiF+whaloF(-1)))^(lambdaw/(1-lambdaw)))=(1-doEHL)*(exp(whaloF)-exp(STEADY_STATE(whaloF)));

//F30 wage setting 1 (EHL)
doEHL*(-exp(FwF) + exp(psiF)/lambdaw*exp(whaloF)^( lambdaw/(lambdaw-1))*exp(lF) + betta*xiw*(exp(wF(+1))/exp(wF))*(exp(pitildewpiF_tp1))^( 1/(1-lambdaw) )*exp(FwF(+1))    )
=(1-doEHL)*(exp(FwF)-exp(STEADY_STATE(FwF)));

//F31 Law of motion of employment or EHL wage setting 2  
doEHL*(-exp(KwF) + ( exp(whaloF)^( lambdaw/(lambdaw-1))*exp(lF))^(1+sigmaL) + betta*xiw*( exp(pitildewpiF_tp1) )^( lambdaw*(1+sigmaL)/(1-lambdaw) )*exp(KwF(+1)))    
=(1-doEHL)*(-exp(lF)+(rho+exp(xF))*exp(lF(-1)));

//F32 wage setting 3 (EHL)
doEHL*(1-(1-xiw)*(AEHL*exp(KwF)/(exp(FwF)*exp(wF)))^( 1/( 1-lambdaw*(1+sigmaL) )) - xiw*(  exp(pitildewpiF)   )^(1/(1-lambdaw)))=(1-doEHL)*(exp(KwF)-exp(STEADY_STATE(KwF)));
//////////////////////////
//Exogenous Variables/////
//////////////////////////

//E1 Composite technology growth
muF=alfa/(1-alfa)*mupsiF+muzF;

//E2 Unit root invest. Tech.
mupsiF=(1-rhomupsi)*STEADY_STATE(mupsiF)+rhomupsi*mupsiF(-1)+sig_mupsi*mupsi_eps/100;

//E3 Unit root neutral Tech.
muzF=(1-rhomuz)*STEADY_STATE(muzF)+rhomuz*muzF(-1)+sig_muz*muz_eps/100;

//E4 Diffusion of composite technology into gov. spending
(1-dolagZBAR)*(exp(nGF)-(exp(nGF(-1))/exp(muF))^(1-thetaG))=dolagZBAR*(exp(nGF)-(exp(nGF(-1)))^(1-thetaG)/exp(muF));

//E5  Diffusion of composite technology into fixed costs of production
nPHIF=nGF;

//E6  Diffusion of composite technology into firm delay costs of bargaining 
nGAMF=nGF;

//E7  Diffusion of composite technology into unemployment benefits
nDF=nGF;

//E8 Diffusion of composite technology into hiring/search costs 
nKAPF=nGF;
////////////////////////////////////////////////////////////////////// 
//Aggregating 'F' and 'R' economies, expressed in percent deviations// 
//////////////////////////////////////////////////////////////////////
//A1 
GDPAGG-STEADY_STATE(GDPF)=GDPR-STEADY_STATE(GDPF)+GDPF-STEADY_STATE(GDPF);
//A2 
piAGG-STEADY_STATE(piF)=piR-STEADY_STATE(piF)+piF-STEADY_STATE(piF);
//A3 
RAGG-STEADY_STATE(RF)=RR-STEADY_STATE(RF)+RF-STEADY_STATE(RF);
//A4 
ukAGG-STEADY_STATE(ukF)=ukR-STEADY_STATE(ukF)+ukF-STEADY_STATE(ukF);
//A5 
lAGG-STEADY_STATE(lF)=lR-STEADY_STATE(lF)+lF-STEADY_STATE(lF);
//A6 
wAGG-STEADY_STATE(wF)=wR-STEADY_STATE(wF)+wF-STEADY_STATE(wF);
//A7 
cAGG-STEADY_STATE(cF)=cR-STEADY_STATE(cF)+cF-STEADY_STATE(cF);
//A8 
iAGG-STEADY_STATE(iF)=iR-STEADY_STATE(iF)+iF-STEADY_STATE(iF);
//A9 
unempAGG-STEADY_STATE(unempF)=unempR-STEADY_STATE(unempF)+unempF-STEADY_STATE(unempF);
//A10 
vTotAGG-STEADY_STATE(vTotF)=vTotR-STEADY_STATE(vTotF)+vTotF-STEADY_STATE(vTotF);
//A11 
fAGG-STEADY_STATE(fF)=fR-STEADY_STATE(fF)+fF-STEADY_STATE(fF);
//A12
pinvestAGG=-mupsiF; 
end; 
/////////////////////////////////
//End Model            //////////
///////////////////////////////// 


% model switches
doNash=0;    %if =0, alt offer barg. ; if =1, Nash bargaining; note: requires doEHL=0 below
doEHL=0;     %if=0, alt offer barg or Nash barg; if=1, EHL labor market
@#define do_given_bets=0  // 1: given values for beta's in sharing rule
%AOB sharing rule coefficients
bet1 = 0.0907;
bet2 = 28.9219;
bet3 = 0.4562;

%Labor market parameters
u=0.055;                 %unemp. rate
rho=0.9;                 %job survival rate
sigm=0.5570;             %matching function share of unemp.
recSHAREpercent=0.5;     %hiring cost para; hiring cost relative to output, in percent
searchSHAREpercent=0.05; %search cost para; hiring cost relative to output, in percent

DSHARE=0.6682;           %unempl. benefits as share of w (replacement ratio)
Q=0.7;                   %vacancy filling rate

deltapercent=0.3022;     %prob. of barg. session determination
M=60;                    %maximum bargaining rounds per quarter, needs to be an even number!!!!
 
if M<2, error('M must be at least equal to 2!');end
if mod(M,2)==1, error('M must be an even number so that the worker makes the last offer!');end 

%prices  
xi=0.5841;         %Calvo prices
pibar=1.00625;     %inflation rate, gross, quarterly
kappaf=0;          %price indexation to past inflation
varkappaf=1;       %price indexation parameter; if kappaf=0 and varkappaf=1 and pibreve=1 -> no indexation at all.
pibreve=1;         %price indexation to value pibreve

lambda=1.4256;     %steady state price markup
nuf=1;             %working capital fraction
tau=1;             %steady state markup shifter

%technology and adjustment cost
alfa=0.2366;       %capital share
deltak=0.025;      %depreciation rate of capital
Spp=13.6684;       %second derivative of invest. adjustment costs
sigmaa=0.0760;     %capacity utilization costs

%growth rates of investment and real GDP
idiffdataPercent=2.9; 
ydiffdataPercent=1.7; 

%Preferences
betta=0.996783170280770; %discount factor households; implies 3% real rate
b=0.8320;                %habit formation in consumption

%Monetary Policy by Taylor rule
rhoR	=	0.8555;   %Interest rate smoothing
rpi     =	1.3552;   %Coefficient on inflation
ry      =	0.0359;   %Coef. on GDP

%Government
etag=0.2;           %Government spending-to-output in ss

%technology diffusion into gov. spending and other parameters of the model
thetaG =0.0136;     %gov spending
thetaPHI=thetaG;    %fixed cost of production
thetaKAP=thetaG;    %hiring cost
thetaGAM=thetaG;	%cost of counteroffer
thetaD=0.7416;      %replacement ratio
dolagZBAR=1;        %if =1, diffusion speed: zbar_t=(zplus_t-1)^thetaj*(zbar_t-1)^(1-thetaj); 
                    %if =0, zbar_t=(zplus_t)^thetaj*(zbar_t-1)^(1-thetaj) for j=G,GAM,KAP,BU,PHI
                    %check cet_steadystate.m for further restrictions on
                    %theta parameters!

% steady state profits as share of output
profy=0;

%EHL parameter
lambdaw=1.2;  %wage markup (EHL)
xiw=0.75;     %Calvo wage stickiness (EHL)
sigmaL=1;     %inv. labor supply elast. (EHL)
kappaw=0;     %wage indexation to past inflation (EHL)  
varkappaw=1;  %wage indexation parameter; if kappaw=0 and varkappaw=1 and pibreve=1 and thetaw=0 -> no indexation at all.
thetaw=0;     %wage indexation to technology

    
% standard deviation of shocks (%)
sig_muz     =	0.1362;       %unit root neutral tech.
sig_mupsi	=	0.1100;       %unit root investment tech.
sig_epsR	=	0.6028;       %monetary policy

sig_epsilon	=	0;            %stationary neutral tech.
sig_upsilon	=	0;            %stationary invest. tech.
sig_g	    =	0;            %gov. spending
sig_taud	=	0;            %price markup
sig_zetac	=	0;            %cons. preference

%AR(1)
rhomupsi	=	0.7279; %unit root investment tech.
rhomuz      =	0;      %unit root neutral tech.

rhoepsilon	=	0;      %stationary neutral tech. 
rhoupsilon	=	0;      %stationary invest. tech.
rhozetac	=	0;      %cons. preference
rhog        =	0;      %gov. spending
rhotaud     =   0;      %price markup
rhosigmam	=	0;      %matching function shock

iota=1;      %This parameter does not apprear in the model anymore. Just ignore it.


shocks; 
var epsR_eps        = 1; //monetary policy
var muz_eps         = 1; //neutral tech.
var mupsi_eps       = 1; //invest. tech.
end;    

steady_state_model;
%% based on cet_steadystate.m, some code is put into cet_steady_helper.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%force the below parameters to be identical to thetaG
thetaPHI=thetaG;%thetaPHI =1; %fixed cost of production
thetaKAP=thetaG;%thetaKAP =1; %hiring cost
thetaGAM=thetaG;%thetaGAM =1; %counteroffer cost
thetaD=thetaG;%thetaD=1;      %unemp. benefits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first, some easy steady states
delta=deltapercent/100;
mu=exp(ydiffdataPercent/100/4);
mupsi=exp(idiffdataPercent/100/4)/mu;
muz=mu/mupsi^(alfa/(1-alfa));
[info,nG,nKAP,nD,nGAM,nPHI] = cet_steady_helper(0,0,0,0,0,  dolagZBAR,mu,thetaG,thetaKAP,thetaD,thetaGAM,thetaPHI);
uk=1;
aofuk=0;
Stilde=0;
Sprimetilde=0;
l=1-u;
R=mu*pibar/betta;
pkprime=1;
Rk=pibar*mu/betta;
%Some more complicated steady states
sigmab=Rk*mupsi*pkprime/pibar-(1-deltak)*pkprime;
aofukprime=sigmab;
pitilde=pibar^kappaf*pibar^(1-kappaf-varkappaf)*pibreve^varkappaf;
mc=1/lambda*(1-betta*xi*(pitilde/pibar)^(lambda/(1-lambda)))/(1-betta*xi*(pitilde/pibar)^(1/(1-lambda))) *((1-xi*(pitilde/pibar)^(1/(1-lambda)))/(1-xi))^(1-lambda);
phalo=((1-xi*(pitilde/pibar)^(1/(1-lambda)))/(1-xi))^(1-lambda)/((1-xi*(pitilde/pibar)^(lambda/(1-lambda)))/(1-xi))^((1-lambda)/lambda);
kl=(alfa*(mupsi*mu)^(1-alfa)*mc/sigmab/tau)^(1/(1-alfa));
%EHL as special case. Obviously, whalo has no meaning otherwise and is set to one!
piw=mu*pibar;
pitildew=(pibar)^kappaw*(pibar)^(1-kappaw-varkappaw)*pibreve^varkappaw*mu^thetaw;
whalo=doEHL*(((1-xiw*(pitildew/piw)^(1/(1-lambdaw)))/(1-xiw))^(1-lambdaw) / ((1-xiw*(pitildew/piw)^(lambdaw/(1-lambdaw)))/(1-xiw))^((1-lambdaw)/lambdaw))+(1-doEHL)*1;
varthet=(1-alfa)*mc/(tau*(mupsi*mu)^alfa*(nuf*R+1-nuf))*kl^alfa;
y=mc/( (phalo^(lambda/(1-lambda))-1)*mc+1-profy)*(kl/(mu*mupsi))^alfa*l*whalo^(lambdaw/(lambdaw-1));
k=kl*l*whalo^(lambdaw/(lambdaw-1));
phi=((k/(l*whalo^(lambdaw/(lambdaw-1)))/(mupsi*mu))^alfa*l*whalo^(lambdaw/(lambdaw-1))-y*phalo^(lambda/(1-lambda)))/nPHI;
i=(1-(1-deltak)/(mu*mupsi))*k;
c=(1-etag)*y-i-(1-doEHL)*recSHAREpercent/100*y-(1-doEHL)*searchSHAREpercent/100*y;
psi=(c-b*c/mu)^(-1)-betta*b*(c*mu-b*c)^(-1);
g=etag*y/nG;
gdp=g*nG+c+i;
K=lambda*psi*y*mc/(1-betta*xi*(pitilde/pibar)^(lambda/(1-lambda)));
F=psi*y/(1-betta*xi*(pitilde/pibar)^(1/(1-lambda)));
%steady states specific to unemployment models, except for real wage (switch)
x=1-rho;
f=x*l/(1-rho*l);
v=x/Q;
sigmam=x*l*(1-rho*l)^(-sigm)*(l*v)^(sigm-1);
kappa=recSHAREpercent/100/x/l*y/nKAP;
s=(searchSHAREpercent/100/((Q^(-1)*x))/l*y)/nKAP;
J=kappa*nKAP+s*nKAP*Q^(-1);
varthetp=varthet/(1-rho*betta);
wp=varthetp-J;
w=doEHL*(varthet)+(1-doEHL)*(wp*(1-rho*betta)); %EHL switch for wage
D=DSHARE*w/nD;
VminusU=(wp-nD*D/(1-rho*betta))/(1-(1-f)*betta*rho)*(1-rho*betta);
U=(nD*D+betta*f*VminusU)/(1-betta);
V=VminusU+U;
A=V-wp;
%AOB sharing rule
alp1=1-delta+(1-delta)^M;
alp2=1-(1-delta)^M;
alp3=alp2*(1-delta)/delta-alp1;
alp4=(1-delta)/(2-delta)*alp2/M+1-alp2;
@#if do_given_bets==0
    bet1=alp2/alp1;
    bet2=alp3/alp1;
    bet3=alp4/alp1;
@#endif
gamma=(bet1*(V-U)-J+bet3*(varthet-nD*D))/(nGAM*bet2);
recruiting=nKAP*kappa*x*l+nKAP*s*(Q^(-1)*x)*l;
%Worker surplus share under Nash sharing
eta=(V-U)/(V-U+J);
%EHL sticky wages
Kw=((l*whalo^(lambdaw/(lambdaw-1)))^(1+sigmaL))/(1-betta*xiw*(pitildew/piw)^(lambdaw*(1+sigmaL)/(1-lambdaw)));
Fw=(psi/lambdaw*(l*whalo^(lambdaw/(lambdaw-1))))/(1-betta*xiw*(pitildew/piw)^(1/(1-lambdaw)));
AEHL=1/(Kw)*(Fw*w)*((1- xiw*(  pitildew/piw  )^(1/(1-lambdaw)))/(1-xiw))^(( 1-lambdaw*(1+sigmaL) ));
%parameters for dynare
tau_SS=log(tau); epsilon_SS=log(1); upsilon_SS=log(1);
zetac_SS=log(1); g_SS=log(g); unemp_SS=log(1-l); vTot_SS=log(v*l);
%steady state variables for ys vector (R economy)
psiR=log(psi); cR=log(c); RR=log(R);  piR=log(pibar);
iR=log(i); pkprimeR=log(pkprime); FR=log(F); KR=log(K);
ukR=log(uk); kbarR=log(k); lR=log(l); xR=log(x); varthetR=log(varthet);
phaloR=log(phalo); yR=log(y); VR=log(V); JR=log(J);   UR=log(U);
GDPR=log(gdp); fR=log(f); vR=log(v); unempR=log(u);
vTotR=log(v*l); QR=log(Q); piwR=log(piw); whaloR=log(whalo);
FwR=log(Fw); KwR=log(Kw); wR=log(w);  wpR=log(wp);
varthetpR=log(varthetp); AAR=log(A);
%steady state variables for ys vector (R economy)
psiF=log(psi); cF=log(c); RF=log(R);  piF=log(pibar);
iF=log(i); pkprimeF=log(pkprime); FF=log(F); KF=log(K);
ukF=log(uk); kbarF=log(k); lF=log(l); xF=log(x); varthetF=log(varthet);
phaloF=log(phalo); yF=log(y); VF=log(V); JF=log(J);   UF=log(U);
GDPF=log(gdp); fF=log(f); vF=log(v); unempF=log(u);  nDF=log(nD);
vTotF=log(v*l); QF=log(Q); piwF=log(piw);nPHIF=log(nPHI);wF=log(w);
whaloF=log(whalo); FwF=log(Fw); KwF=log(Kw); muzF=log(muz);
mupsiF=log(mupsi); muF=log(mu); nGF=log(nG); nGAMF=log(nGAM); nKAPF=log(nKAP);
AAF=log(A);wpF=log(wp); varthetpF=log(varthetp);
%AGG econ
GDPAGG=log(gdp); piAGG=log(pibar); RAGG=log(R); ukAGG=log(uk); lAGG=log(l);
wAGG=wF; cAGG=log(c); iAGG=log(i); unempAGG=log(u); vTotAGG=log(v*l);
fAGG=log(f); pinvestAGG=-log(mupsi);

info = cet_steady_helper(1,Kw,J,V,U); % THIS CHECKS WHETHER STEADY-STATES ARE NONZERO
end;

resid;
steady;
check;


@#ifdef ML
%% FREQUENTIST ESTIMATION
estimated_params;
xi;
lambda             ,     , 1.001, Inf;
rhoR               , 0.85,  -Inf, Inf;
rpi                , 1.45,  -Inf, Inf;
ry                 , 0.01,  -Inf, Inf;
b                  ,     ,  -Inf, Inf;
sigmaa             ,     ,  -Inf, Inf;
Spp                ,     ,  -Inf, Inf;
alfa               ,     ,  -Inf, 0.6;
thetaG             , 0.10,     0,   1;
deltapercent       ,     ,  -Inf, Inf;
DSHARE             ,     ,  -Inf, Inf;
recSHAREpercent    ,     ,  -Inf, Inf;
searchSHAREpercent ,     ,  -Inf, Inf;
sigm               ,     ,  -Inf, Inf;
sig_epsR           ,     ,  -Inf, Inf;
sig_muz            ,     ,  -Inf, Inf;
sig_mupsi          ,     ,  -Inf, Inf;
rhomupsi           ,     ,  -Inf, Inf;
end;

xi = 0.7363; 
lambda = 1.4082;
rhoR = 0.8463;
rpi = 1.4566;
ry = 0.0407;
b = 0.8196;
sigmaa = 0.1544;
Spp = 15.1062;
alfa = 0.2404;
thetaG = 0.0502;
deltapercent = 0.1825;
DSHARE = 0.4596;
recSHAREpercent = 0.4849;
searchSHAREpercent = 0.0654;
sigm = 0.5403;
sig_epsR = 0.6163;
sig_muz = 0.1496;
sig_mupsi = 0.1150;
rhomupsi = 0.7314;

estimated_params_init(use_calibration);
end;

@#else

%% BAYESIAN ESTIMATION
estimated_params;
xi                 ,     ,      , , beta_pdf  , 0.66, 0.100;
lambda             ,     , 1.001, , gamma_pdf , 1.20, 0.050;
rhoR               , 0.85,      , , beta_pdf  , 0.70, 0.150;
rpi                , 1.45,      , , gamma_pdf , 1.70, 0.150;
ry                 , 0.01,      , , gamma_pdf , 0.10, 0.050;
b                  ,     ,      , , beta_pdf  , 0.50, 0.150;
sigmaa             ,     ,      , , gamma_pdf , 0.50, 0.300;
Spp                ,     ,      , , gamma_pdf , 8.00, 2.000;
alfa               ,     ,      , , beta_pdf  , 0.33, 0.025;
thetaG             , 0.10,      , , beta_pdf  , 0.50, 0.200;
deltapercent       ,     ,      , , gamma_pdf , 0.50, 0.400;
DSHARE             ,     ,      , , beta_pdf  , 0.40, 0.100;
recSHAREpercent    ,     ,      , , gamma_pdf , 1.00, 0.300;
searchSHAREpercent ,     ,      , , gamma_pdf , 0.10, 0.070;
sigm               ,     ,      , , beta_pdf  , 0.50, 0.100;
sig_epsR           ,     ,      , , gamma_pdf , 0.65, 0.050;
sig_muz            ,     ,      , , gamma_pdf , 0.10, 0.050;
sig_mupsi          ,     ,      , , gamma_pdf , 0.10, 0.050;
rhomupsi           ,     ,      , , beta_pdf  , 0.75, 0.100;
end;
@#endif

varobs GDPAGG piAGG RAGG ukAGG lAGG wAGG cAGG iAGG pinvestAGG unempAGG vTotAGG fAGG;