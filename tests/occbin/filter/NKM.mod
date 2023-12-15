//Tests Occbin estimation with IVF and PKF with 1 constraints

// this file implements the model in:
// Atkinson, T., A. W. Richter, and N. A. Throckmorton (2019). 
// The zero lower bound andestimation accuracy.Journal of Monetary Economics
// original codes provided by Alexander Richter
// adapted for dynare implementation
// ------------------------- Settings -----------------------------------------//

@#ifndef small_model
   @#define small_model = 0
@#endif

// ----------------- Defintions -----------------------------------------//
var        
    c       $c$         (long_name='Consumption')
    n                   (long_name='Labor')
    y       $y$         (long_name='Output')
    yf                  (long_name='Final goods')
    yg                  (long_name='Output growth gap')
    w                   (long_name='Real wage rate')
    wf                  (long_name='Flexible real wage')
    pigap               (long_name='Inflation rate -> pi(t)/pibar = pigap')
    inom    ${i^{nom}}$ (long_name='Nominal interest rate')
    inomnot             (long_name='Notional interest rate')
    mc                  (long_name='Real marginal cost')
    lam  ${\lambda}$    (long_name='Inverse marginal utility of wealth')
    g                   (long_name='Growth shock')
    s                   (long_name='Risk premium shock')
    mp                  (long_name='Monetary policy shock')
    pi   ${\pi}$        (long_name='Observed inflation')
    @#if !(small_model)
        x               (long_name='Investment')
        k               (long_name='Capital')
        u               (long_name='Utilization cost')
        ups             (long_name='Utilization choice')
        wg              (long_name='Real wage growth gap')    
        xg              (long_name='Investment growth')
        rk              (long_name='Real rental rate')
        q               (long_name='Tobins q')
    @#endif
;
varexo          
    epsg    ${\varepsilon_g}$   (long_name='Productivity growth shock')
    epsi                        (long_name='Notional interest rate shock')
    epss                        (long_name='Risk premium shock')
;        
parameters
    // Calibrated Parameters    
    beta    $\beta$    (long_name='Discount factor')
    chi                 (long_name='Labor disutility scale')
    thetap              (long_name='Elasticity of subs. between intermediate goods')
    thetaw              (long_name='Elasticity of subs. between labor types')
    nbar                (long_name='Steady state labor')
    eta                 (long_name='Inverse frish elasticity of labor supply')
    delta               (long_name='Depreciation')
    alpha               (long_name='Capital share')
    gbar                (long_name='Mean growth rate')
    pibar               (long_name='Inflation target')
    inombar             (long_name='Steady gross nom interest rate')
    inomlb              (long_name='Effective lower bound on gross nominal interest rate')
    sbar                (long_name='Average risk premium')
    // Parameters for DGP and Estimated parameters
    varphip             (long_name='Rotemberg price adjustment cost')
    varphiw             (long_name='Rotemberg wage adjustment cost')
    h          	        (long_name='Habit persistence')
    rhos                (long_name='Persistence')
    rhoi    	        (long_name='Persistence')
    sigz                (long_name='Standard deviation technology')
    sigs                (long_name='Standard deviation risk premia')
    sigi                (long_name='Standard deviation mon pol')
    phipi               (long_name='Inflation responsiveness')
    phiy                (long_name='Output responsiveness')
    nu                  (long_name='Investment adjustment cost')
    sigups              (long_name='Utilization')
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
    yg = g*y/(gbar*y(-1));
       
    [name = 'Notional Interest Rate (9)']
    inomnot = inomnot(-1)^rhoi*(inombar*pigap^phipi*yg^phiy)^(1-rhoi)*exp(mp);    
    
    [name = 'Nominal Interest Rate (10)', bind='zlb']
    inom = inomlb;
    [name = 'Nominal Interest Rate (10)', relax='zlb']
    inom = inomnot;
                  
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
options_.TeX=1;
occbin_constraints;
name 'zlb'; bind inom <=  inomlb; relax inom > inomlb;
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
    options_.occbin.smoother.debug=1;
    occbin_setup(filter_use_relaxation,likelihood_piecewise_kalman_filter);
    // use PKF  
    estimation(
            datafile=dataobsfile2, mode_file=NKM_mh_mode_saved,
            mode_compute=0, nobs=120, first_obs=1,
            mh_replic=0, plot_priors=0, smoother,
            consider_all_endogenous,heteroskedastic_filter,filter_step_ahead=[1],smoothed_state_uncertainty);
    
    // plot regimes
    occbin.plot_regimes(oo_.occbin.smoother.regime_history,M_,options_)

    // forecast starting from period 42, zero shocks (default)
    smoother2histval(period=42);
    [forecast, error_flag] = occbin.forecast(options_,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,8);
    // forecast with stochastic shocks
    options_.occbin.forecast.qmc=true;
    options_.occbin.forecast.replic=127;
    [forecast1, error_flag] = occbin.forecast(options_,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,8);

    figure('Name','OccBin: Forecasts')
    subplot(2,1,1)
    plot(1:8,forecast.piecewise.Mean.inom,'b-',1:8,forecast1.piecewise.Mean.inom,'r--')
    subplot(2,1,2)
    plot(1:8,forecast.piecewise.Mean.y,'b-',1:8,forecast1.piecewise.Mean.y,'r--')

    // GIRF given states in 42 and shocks in 43
    t0=42;
    options_.occbin.irf.exo_names=M_.exo_names;
    options_.occbin.irf.t0=t0;
    oo_.occbin.irfs = occbin.irf(M_,oo_,options_);

    var_list_ = {'c','n','y','pigap','inom','inomnot'};

    // if you want to scale ...
    // options_occbin_.plot_irf.endo_scaling_factor = vars_irf(:,3);
    options_.occbin.plot_irf.simulname = ['t0_' int2str(t0)];
    options_.irf=40;
    occbin.plot_irfs(M_,oo_.occbin.irfs,options_,var_list_);
    var_list_={};
    options_.occbin.plot_irf.simulname = ['t0_' int2str(t0) '_full'];
    occbin.plot_irfs(M_,oo_.occbin.irfs,options_,var_list_);
    oo0=oo_;
    // use smoother_redux
    estimation(
            datafile=dataobsfile2, mode_file=NKM_mh_mode_saved,
            mode_compute=0, nobs=120, first_obs=1,
            mh_replic=0, plot_priors=0, smoother, smoother_redux,
            consider_all_endogenous,heteroskedastic_filter,filter_step_ahead=[1],smoothed_state_uncertainty);

    // check consistency of smoother_redux
    for k=1:M_.endo_nbr, 
        mer(k)=max(abs(oo_.SmoothedVariables.(M_.endo_names{k})-oo0.SmoothedVariables.(M_.endo_names{k}))); 
    end
    if max(mer)>1.e-10
        error('smoother redux does not recover full smoother results!')
    else
        disp('smoother redux successfully recovers full smoother results!')
    end

    // use inversion filter (note that IF provides smoother together with likelihood)
    occbin_setup(likelihood_inversion_filter,smoother_inversion_filter);
            
    estimation(
            datafile=dataobsfile2, mode_file=NKM_mh_mode_saved,
            mode_compute=0, nobs=120, first_obs=1,
            mh_replic=0, plot_priors=0, smoother,
            consider_all_endogenous,heteroskedastic_filter,filter_step_ahead=[1],smoothed_state_uncertainty);
            
    // show initial condition effect of IF
    figure('Name','OccBin: Smoothed shocks')
    subplot(221)
    plot([oo0.SmoothedShocks.epsg oo_.SmoothedShocks.epsg]), title('epsg')
    subplot(222)
    plot([oo0.SmoothedShocks.epsi oo_.SmoothedShocks.epsi]), title('epsi')
    subplot(223)
    plot([oo0.SmoothedShocks.epss oo_.SmoothedShocks.epss]), title('epss')
    legend('PKF','IF')
    figure('Name','OccBin: Smoothed Variables')
    subplot(221)
    plot([oo0.SmoothedVariables.inom oo_.SmoothedVariables.inom]), title('inom')
    subplot(222)
    plot([oo0.SmoothedVariables.yg oo_.SmoothedVariables.yg]), title('yg')
    subplot(223)
    plot([oo0.SmoothedVariables.pi oo_.SmoothedVariables.pi]), title('pi')
    legend('PKF','IF')
    occbin_write_regimes(smoother);

write_latex_dynamic_model;
collect_latex_files;
[status, cmdout]=system(['pdflatex -halt-on-error -interaction=nonstopmode ' M_.fname '_TeX_binder.tex']);
if status
    cmdout
    error('TeX-File did not compile.')
end
