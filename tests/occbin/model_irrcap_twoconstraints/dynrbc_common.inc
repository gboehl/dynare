steady_state_model;
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;


k = log(kss);
c = log(css);
i = log(iss);
lambdak = 0;
a=0;
end;

BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;
GAMMAC=2;
RHOA = 0.9;
PHI = 0.975;
PSI = 5;        % adjustment cost for capital if investment is negative

shocks;
  var erra; stderr 0.015;
end;

steady;

// run occbin simulations
% Option=1: impulse responses
% Option=2: random simulation

@#for option_val in [1, 2]

option=@{option_val};

%%%%%%%%%%%%%%%% Inputs stop here %%%%%%%%%%%%%%%%%%%%%

if option==1
   shocks(surprise,overwrite);
   var erra;
   periods 1:9, 10, 50, 90, 130, 131:169;
   values -0.0001, -0.01,-0.02, 0.01, 0.02, 0;
   end;


elseif option==2
   nperiods = 40;
   randn('seed',1);
   shockssequence = 1*randn(nperiods,1)*0.02 ;

   shocks(surprise,overwrite);
   var erra;
   periods 1:40;
   values (shockssequence);
   end;

end

% set inputs

occbin_setup(smoother_debug);
occbin_solver(simul_debug,simul_periodic_solution,simul_periods=200,simul_maxit=300,simul_curb_retrench,simul_check_ahead_periods=300);
occbin_write_regimes(filename='test',periods=[1:40]);

%% Modify to plot IRFs
titlelist = char('c','lambdak','k','i','a');
percent = 'Percent';
value = 'value';
ylabels = char(percent,value,percent,percent,percent);
figtitle = 'Simulated variables';
legendlist = cellstr(char('Piecewise Linear','Linear'));

oo_= occbin.unpack_simulations(M_,oo_,options_);

line1=100*[oo_.occbin.endo_piecewise.c-oo_.occbin.endo_ss.c,oo_.occbin.endo_piecewise.lambdak/100,oo_.occbin.endo_piecewise.k-oo_.occbin.endo_ss.k,oo_.occbin.endo_piecewise.i-oo_.occbin.endo_ss.i,oo_.occbin.endo_piecewise.a-oo_.occbin.endo_ss.a];
line2=100*[oo_.occbin.endo_linear.c-oo_.occbin.endo_ss.c,oo_.occbin.endo_linear.lambdak/100,oo_.occbin.endo_linear.k-oo_.occbin.endo_ss.k,oo_.occbin.endo_linear.i-oo_.occbin.endo_ss.i,oo_.occbin.endo_linear.a-oo_.occbin.endo_ss.a];

occbin.make_chart(titlelist,legendlist,figtitle,ylabels,cat(3,line1,line2));

occbin_graph(noconstant) c erra lambdak k i a k;

@#if option_val==1
    verbatim;
    c=oo_.occbin.endo_piecewise.c;
    end;
    save('datasim.mat','c');
    varobs c;

    occbin_solver(simul_periods=200,simul_maxit=200,simul_curb_retrench,simul_check_ahead_periods=200);
    occbin_setup(smoother_periods=200,smoother_maxit=200,smoother_curb_retrench,smoother_check_ahead_periods=200);
    calib_smoother(datafile=datasim);
    oo0=oo_;
    occbin_setup(filter_use_relaxation);
    calib_smoother(datafile=datasim);
    oo_= occbin.unpack_simulations(M_,oo_,options_);

    titlelist = char('c','lambdak','k','i','a','erra');
    percent = 'Percent';
    value = 'value';
    ylabels = char(percent,value,percent,percent,percent, value);
    figtitle = 'Smoothed variables (piecewise)';
    legendlist = cellstr(char('Simulated','Piecewise smoother'));

    shock_vector=[oo_.occbin.simul.shocks_sequence./100; zeros(length(oo_.occbin.endo_piecewise.c)-size(oo_.occbin.simul.shocks_sequence,1),size(oo_.occbin.simul.shocks_sequence,2))];
    line1=100*[oo_.occbin.endo_piecewise.c-oo_.occbin.endo_ss.c,oo_.occbin.endo_piecewise.lambdak/100,oo_.occbin.endo_piecewise.k-oo_.occbin.endo_ss.k,oo_.occbin.endo_piecewise.i-oo_.occbin.endo_ss.i,oo_.occbin.endo_piecewise.a-oo_.occbin.endo_ss.a, shock_vector];
    line2=100*[oo_.occbin.smoother.SmoothedVariables.c-oo_.occbin.endo_ss.c,oo_.occbin.smoother.SmoothedVariables.lambdak/100,oo_.occbin.smoother.SmoothedVariables.k-oo_.occbin.endo_ss.k,oo_.occbin.smoother.SmoothedVariables.i-oo_.occbin.endo_ss.i,oo_.occbin.smoother.SmoothedVariables.a-oo_.occbin.endo_ss.a, oo_.occbin.smoother.SmoothedShocks.erra/100];
    occbin.make_chart(titlelist,legendlist,figtitle,ylabels,cat(3,line1,line2));

    figtitle = 'Smoothed variables (linear)';
    legendlist = cellstr(char('Simulated','Linear smoother'));

    line1=100*[oo_.occbin.endo_piecewise.c-oo_.occbin.endo_ss.c,oo_.occbin.endo_piecewise.lambdak/100,oo_.occbin.endo_piecewise.k-oo_.occbin.endo_ss.k,oo_.occbin.endo_piecewise.i-oo_.occbin.endo_ss.i,oo_.occbin.endo_piecewise.a-oo_.occbin.endo_ss.a, shock_vector];
    line2=100*[oo_.occbin.linear_smoother.SmoothedVariables.c-oo_.occbin.endo_ss.c,oo_.occbin.linear_smoother.SmoothedVariables.lambdak/100,oo_.occbin.linear_smoother.SmoothedVariables.k-oo_.occbin.endo_ss.k,oo_.occbin.linear_smoother.SmoothedVariables.i-oo_.occbin.endo_ss.i,oo_.occbin.linear_smoother.SmoothedVariables.a-oo_.occbin.endo_ss.a, oo_.occbin.linear_smoother.SmoothedShocks.erra/100];
    occbin.make_chart(titlelist,legendlist,figtitle,ylabels,cat(3,line1,line2));

    occbin.plot_regimes(oo_.occbin.smoother.regime_history,M_,options_)
@#endif
@#endfor
