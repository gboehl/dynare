//Tests Occbin estimation with IVF and PKF with 1 constraints and redundant shocks

// this file implements the model in:
// Atkinson, T., A. W. Richter, and N. A. Throckmorton (2019). 
// The zero lower bound andestimation accuracy.Journal of Monetary Economics
// original codes provided by Alexander Richter
// adapted for dynare implementation
// ------------------------- Settings -----------------------------------------//

@#ifndef small_model
   @#define small_model = 0
@#endif

// if ~exist('run_ivf','var')
run_ivf=0;
// end        

// ----------------- Defintions -----------------------------------------//
var        
    c          //1  Consumption 
    n          //2  Labor
    y          //5  Output
    yf         //6  Final goods       
    yg         //11 Output growth gap
    w          //12 Real wage rate
    wf         //13 Flexible real wage
    pigap      //15 Inflation rate -> pi(t)/pibar = pigap
    inom       //16 Nominal interest rate
    inomnot    //17 Notional interest rate
    mc         //19 Real marginal cost
    lam        //20 Inverse marginal utility of wealth  
    g          //21 Growth shock       
    s          //22 Risk premium shock
    mp         //23 Monetary policy shock    
    pi      //24 Observed inflation    
    @#if !(small_model)
        x          //3  Investment
        k          //4  Capital    
        u          //7  Utilization cost
        ups        //8  Utilization choice
        wg         //9  Real wage growth gap    
        xg         //10 Investment growth
        rk         //14 Real rental rate
        q          //18 Tobins q   
    @#endif
;
varexo          
    junk1
    epsg       // Productivity growth shock
    epsi       // Notional interest rate shock
    epss       // Risk premium shock
    junk2
;        
parameters
    // Calibrated Parameters    
    beta        // Discount factor
    chi         // Labor disutility scale
    thetap      // Elasticity of subs. between intermediate goods
    thetaw      // Elasticity of subs. between labor types
    nbar        // Steady state labor
    eta         // Inverse frish elasticity of labor supply
    delta       // Depreciation
    alpha       // Capital share
    gbar        // Mean growth rate
    pibar       // Inflation target
    inombar     // Steady gross nom interest rate    
    inomlb      // Effective lower bound on gross nominal interest rate    
    sbar        // Average risk premium
    // Parameters for DGP and Estimated parameters
    varphip     // Rotemberg price adjustment cost
    varphiw     // Rotemberg wage adjustment cost
    h          	// Habit persistence
    rhos        // Persistence
    rhoi    	// Persistence
    sigz        // Standard deviation technology
    sigs        // Standard deviation risk premia
    sigi        // Standard deviation mon pol
    phipi       // Inflation responsiveness
    phiy        // Output responsiveness
    nu          // Investment adjustment cost
    sigups      // Utilization    
        
    // Switching parameters
    zlb 
 ;


// ---------------- Calibration -----------------------------------------//   

beta     = 0.9949;    // Discount factor
thetap   = 6;         // Elasticity of subs. between intermediate goods
thetaw   = 6;         // Elasticity of subs. between labor types
nbar     = 1/3;       // Steady state labor
eta      = 1/3;       // Inverse frish elasticity of labor supply
delta    = 0.025;     // Depreciation
alpha    = 0.35;      // Capital share
gbar     = 1.0034;    // Mean growth rate
pibar    = 1.0053;    // Inflation target
sbar     = 1.0058;    // Average risk premium
rkbar    = 1/beta;
varphiw  = 100;       // Rotemberg wage adjustment cost
nu       = 4;         // Investment adjustment cost
sigups   = 5;         // Utilization
varphip  = 100;       // Rotemberg price adjustment cost
h        = 0.80;      // Habit persistence
rhos     = 0.80;      // Persistence
rhoi     = 0.80; 	  // Persistence
sigz     = 0.005;     // Standard deviation
sigs     = 0.005;     // Standard deviation
sigi     = 0.002;      // Standard deviation
phipi    = 2.0;       // Inflation responsiveness
phiy     = 0.5;       // Output responsiveness
zlb      = 0 ;         // ZLB dummy
inomlb   = 1 ;         // Inom LB  
        
// ---------------- Model -----------------------------------------------//               
model;
    
    @#if !(small_model)
        [name = 'HH FOC utilization (1)']
        rk = steady_state(rk)*exp(sigups*(ups-1));

        [name = 'Utilization definition (3)']
        u = steady_state(rk)*(exp(sigups*(ups-1))-1)/sigups;

        [name = 'Firm FOC capital (4)']
        rk = mc*alpha*g*yf/(ups*k(-1));

        [name = 'HH FOC capital (17)']
        q = beta*(lam/lam(+1))*(rk(+1)*ups(+1)-u(+1)+(1-delta)*q(+1))/g(+1);

        [name = 'HH FOC investment (18)']
        1 = q*(1-nu*(xg-1)^2/2-nu*(xg-1)*xg)+beta*nu*gbar*q(+1)*(lam/lam(+1))*xg(+1)^2*(xg(+1)-1)/g(+1);

        [name = 'Wage Phillips Curve (20)']
        varphiw*(wg-1)*wg = ((1-thetaw)*w+thetaw*wf)*n/yf + beta*varphiw*(lam/lam(+1))*(wg(+1)-1)*wg(+1)*(yf(+1)/yf)   ;

        [name = 'Real wage growth gap (6)']
        wg = pigap*g*w/(gbar*w(-1));
        
        [name = 'Law of motion for capital (15)']
        k = (1-delta)*(k(-1)/g)+x*(1-nu*(xg-1)^2/2);                   

        [name = 'Investment growth gap (14)']
        xg = g*x/(gbar*x(-1));    
    @#endif 
    
            
    @#if small_model         
        [name = 'Production function (2)']
        yf = n;   
    
        [name = 'Firm FOC labor (5)']            
        w = mc*yf/n;
        
        [name = 'Output definition (7)']
        y = (1-varphip*(pigap-1)^2/2)*yf;
        
        [name = 'ARC (13)']
        c = y;
        
        [name = 'Household labpur supply equals flex wage']
        w = wf;
    @#else     
        [name = 'Production function (2)']
        yf = (ups*k(-1)/g)^alpha*n^(1-alpha);
    
        [name = 'Firm FOC labor (5)']            
        w = (1-alpha)*mc*yf/n;
    
       	[name = 'Output definition (7)']
        y = (1-varphip*(pigap-1)^2/2-varphiw*(wg-1)^2/2)*yf - u*k(-1)/g;
      
        [name = 'ARC (13)']
        x = y-c;
    @#endif 
    
    [name = 'Output growth gap (8)']
    yg = g*y/(gbar*y(-1)) + junk1 + junk2;
       
    [name = 'Notional Interest Rate (9)']
    inomnot = inomnot(-1)^rhoi*(inombar*pigap^phipi*yg^phiy)^(1-rhoi)*exp(mp);    
    
    [name = 'Nominal Interest Rate (10)']
    inom = (inomnot*(1-zlb)+zlb*inomlb);   
                  
    [name = 'Inverse MUC (11)']
    lam = c-h*c(-1)/g;
    
    [name = 'Flexible real wage definition (12)']
    wf = chi*n^eta*lam;
    
    [name = 'HH FOC bond (16)']
    1 = beta*(lam/lam(+1))*s*inom/(g(+1)*pibar*pigap(+1));
       
    [name = 'Price Phillips Curve (19)']
    varphip*(pigap-1)*pigap = 1-thetap+thetap*mc+beta*varphip*(lam/lam(+1))*(pigap(+1)-1)*pigap(+1)*(yf(+1)/yf);
     
    [name = 'Stochastic productivity growth (21)']
    g = gbar+ sigz*epsg;
    
    [name = 'Risk premium shock (22)']
    s = (1-rhos)*sbar + rhos*s(-1)+sigs*epss     ;      
    
    [name = 'Notional interest rate shock (23)']
    mp = sigi*epsi;
    
    [name = 'Observed inflation (24)']
    pi = pigap*pibar;
     
end;

occbin_constraints;
name 'zlb'; bind inom+inombar <=  inomlb; relax inom+inombar > inomlb;
end;

// ---------------- Steady state -----------------------------------------//        
steady_state_model;
    mp = 0;
    xg = 1;
    s = sbar;
    n = nbar;
    ups =1;
    u =0;
    q =1;
    g = gbar;
    yg = 1;
    pigap = 1;
    wg =1;
    pi = pigap*pibar;
    // FOC bond 
    inom = gbar*pibar/(beta*s);
    inomnot = inom;
    inombar = inom;
    // Firm pricing
    mc = (thetap-1)/thetap;
    // FOC capital
    rk = gbar/beta+delta-1;
    // Marginal cost definition
    @#if small_model
        alpha = 0;
        thetaw=1;
    @#endif
    w = (mc*(1-alpha)^(1-alpha)*alpha^alpha/rk^alpha)^(1/(1-alpha));
    @#if small_model
        wf = w;
    @#else
        wf = (thetaw-1)*w/thetaw;
    @#endif
    // Consolidated FOC firm
    k = w*n*gbar*alpha/(rk*(1-alpha));
    // Law of motion for capital
    x = (1-(1-delta)/gbar)*k;
    // Production function
    yf = (k/gbar)^alpha*n^(1-alpha);
    // Real GDP
    y = yf;
    // Aggregate resouce constraint
    c = y-x;
    // FOC labor
    lam = (1-h/gbar)*c;
    chi = wf/(n^eta*lam);
    
    // try log observables
end;

// ---------------- Checks -----------------------------------------//        
steady;
check;

// ---------------- Simulation -----------------------------------------//        
shocks;
    var epsi   =  1;
    var epss   =  1;
    var epsg   =  1;
end;

steady;
check;  
        
// ---------------- Estimation -----------------------------------------//        
        
varobs yg inom pi;
    estimated_params;
        // PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
        // PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
        varphip,,0,inf,NORMAL_PDF,100,25;
        phipi,,,,NORMAL_PDF,2,0.25;
        phiy,,0,inf,NORMAL_PDF,0.5,0.25;
        h,,,,BETA_PDF,0.8,0.1;
        rhos,,,,BETA_PDF,0.8,0.1;
        rhoi,,,,BETA_PDF,0.8,0.1;
        sigz,,,,INV_GAMMA_PDF,0.005,0.005;
        sigs,,,,INV_GAMMA_PDF,0.005,0.005;
        sigi,,,,INV_GAMMA_PDF,0.002,0.002;
    end;    
    
 
//     dataloading_jme_beta(1,'sims.txt',30);
    load('dataobsfile','inom')
    // check if inom is at lb and remove data + associated shock
    verbatim;
    inom(inom==1)=NaN;
    end;
    inx = strmatch('epsi',M_.exo_names);
    if any(isnan(inom))
        M_.heteroskedastic_shocks.Qscale_orig.periods=find(isnan(inom));
        M_.heteroskedastic_shocks.Qscale_orig.exo_id=inx;
        M_.heteroskedastic_shocks.Qscale_orig.scale=0;
    else
        options_.heteroskedastic_filter=false;
    end
            
    copyfile dataobsfile.mat dataobsfile2.mat
    save dataobsfile2 inom -append
    // -----------------Occbin ----------------------------------------------//   
    options_.occbin.filter.use_relaxation=true;
    // use PKF  
    estimation(
            datafile=dataobsfile2, mode_file=NKM_mh_mode_saved,
            mode_compute=0, nobs=120, first_obs=1,
            mh_replic=0, plot_priors=0, smoother,
            graph_format=(fig), nodisplay,consider_all_endogenous,heteroskedastic_filter);
    
    oo0=oo_;
    
    // use inversion filter (note that IF provides smoother together with likelihood)
    options_.occbin.likelihood.inversion_filter  = 1;
    options_.occbin.smoother.inversion_filter  = 1;
            
    estimation(
            datafile=dataobsfile2, mode_file=NKM_mh_mode_saved,
            mode_compute=0, nobs=120, first_obs=1,
            mh_replic=0, plot_priors=0, smoother,
            graph_format=(fig), nodisplay, consider_all_endogenous,heteroskedastic_filter);
            
    // show initial condition effect of IF
    figure,
    subplot(221)
    plot([oo0.SmoothedShocks.epsg oo_.SmoothedShocks.epsg]), title('epsg')
    subplot(222)
    plot([oo0.SmoothedShocks.epsi oo_.SmoothedShocks.epsi]), title('epsi')
    subplot(223)
    plot([oo0.SmoothedShocks.epss oo_.SmoothedShocks.epss]), title('epss')
    legend('PKF','IF')
