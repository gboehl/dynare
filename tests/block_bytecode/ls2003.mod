/*
  The model in this file is used to test all combinations of algorithms with
  block and bytecode options.

  It is designed is such a way that it includes the 8 possible types of blocks:
  – Solve {forward, backward, two boundaries} {simple, complete}
  – Evaluate {forward, backward}

  All the “Solve” blocks are also present in both linear and nonlinear forms
  (since the codepaths are typically different depending on the linearity of
  the block).

  Note that there is no such thing as a nonlinear “Evaluate” block, since the
  endogenous variables of the block always enter linearly (on the LHS).
*/

var y y_s R pie dq pie_s de A;
var eb sfc1 sfc2 sbc1 sbc2 sbc3 stbs sfs sbs;
var nsfc1 nsfc2 nsbc1 nsbc2 nsbc3 nstbs nsfs nsbs;
var nstbc_c nstbc_y nstbc_k nstbc_h;

varexo e_R e_q e_ys e_pies e_A;
varexo e_eb e_sfc1 e_sfc2 e_sbc e_stbs e_sbs e_sfs;
varexo e_nsfc1 e_nsfc2 e_nsbc e_nstbs e_nsbs e_nsfs;
varexo e_nstbc;

parameters psi1 psi2 psi3 rho_R tau alpha rr k rho_q rho_A rho_ys rho_pies;
parameters nstbc_beta nstbc_alpha nstbc_delta nstbc_theta nstbc_psi;

psi1 = 1.54;
psi2 = 0.25;
psi3 = 0.25;
rho_R = 0.5;
alpha = 0.3;
rr = 2.51;
k = 0.5;
tau = 0.5;
rho_q = 0.4;
rho_A = 0.2;
rho_ys = 0.9;
rho_pies = 0.7;

nstbc_alpha = 0.36;
nstbc_beta  = 0.99;
nstbc_delta = 0.025;
nstbc_psi   = 0;
nstbc_theta = 2.95;


@#if !block && !bytecode && !use_dll
model;
@#elseif block && !bytecode && !use_dll
model(block, cutoff=0);
@#elseif !block && bytecode
model(bytecode);
@#elseif block && bytecode
model(block, bytecode, cutoff=0);
@#elseif !block && use_dll
model(use_dll);
@#else
model(block, use_dll, cutoff=0);
@#endif

  // Block of type “Solve two boundaries complete” (linear)
  y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
  pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
  pie = de+(1-alpha)*dq+pie_s;
  R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;

  // Block of type “Evaluate forward”
  dq = rho_q*dq(-1)+e_q;
  y_s = rho_ys*y_s(-1)+e_ys;
  pie_s = rho_pies*pie_s(-1)+e_pies;
  A = rho_A*A(-1)+e_A;

  // Block of type “Solve forward complete” (linear)
  sfc1 = 0.2*sfc2+0.5*sfc1(-1)+e_sfc1;
  sfc2 = 0.1*sfc1+0.5*sfc2(-1)+e_sfc2;

  // Block of type “Solve backward complete” (linear)
  sbc1 = sbc1(1)-0.5*(sbc3-sbc2(1))+e_sbc;
  sbc2 = 0.1*sbc1 + 0.9*sbc2(1);
  sbc3 = 1.5*sbc2+0.5*sbc1;

  // Block of type “Evaluate backward”, see #1727
  eb = 0.9*eb(+1) + e_eb;

  // Block of type “Solve two boundaries simple” (linear)
  stbs - 0.8*stbs(-1) = 0.9*(stbs(+1)-0.8*stbs) + e_stbs;

  // Block of type “Solve backward simple” (linear)
  // NB: The LHS is deliberately not a single variable, otherwise is would be classified as “Evaluate backward”
  sbs - 0.8*sbs(1) = e_sbs;

  // Block of type “Solve forward simple” (linear)
  // NB: The LHS is deliberately not a single variable, otherwise is would be classified as “Evaluate forward”
  sfs - 0.8*sfs(-1) = e_sfs;


  // Block of type “Solve forward complete” (nonlinear)
  log(nsfc1) = 0.2*log(nsfc2)+0.5*log(nsfc1(-1))+e_nsfc1;
  log(nsfc2) = 0.1*log(nsfc1)+0.5*log(nsfc2(-1))+e_nsfc2;

  // Block of type “Solve backward complete” (nonlinear)
  log(nsbc1) = log(nsbc1(1))-0.5*(log(nsbc3)-log(nsbc2(1)))+e_nsbc;
  log(nsbc2) = 0.1*log(nsbc1) + 0.9*log(nsbc2(1));
  log(nsbc3) = 1.5*log(nsbc2)+0.5*log(nsbc1);

  // Block of type “Solve two boundaries simple” (nonlinear)
  log(nstbs) - 0.8*log(nstbs(-1)) = 0.9*(log(nstbs(+1))-0.8*log(nstbs)) + e_nstbs;

  // Block of type “Solve backward simple” (nonlinear)
  log(nsbs) = 0.8*log(nsbs(1)) + e_nsbs;

  // Block of type “Solve forward simple” (nonlinear)
  log(nsfs) = 0.8*log(nsfs(-1)) + e_nsfs;

  // Block of type “Solve two boundaries complete” (nonlinear)
  // NB: This is a variation of example1.mod
  nstbc_c*nstbc_theta*nstbc_h^(1+nstbc_psi)=(1-nstbc_alpha)*nstbc_y;
  nstbc_k = nstbc_beta*(((exp(e_nstbc)*nstbc_c)/(exp(e_nstbc(+1))*nstbc_c(+1)))*(exp(e_nstbc(+1))*nstbc_alpha*nstbc_y(+1)+(1-nstbc_delta)*nstbc_k));
  nstbc_y = (nstbc_k(-1)^nstbc_alpha)*(nstbc_h^(1-nstbc_alpha));
  nstbc_k = exp(e_nstbc)*(nstbc_y-nstbc_c)+(1-nstbc_delta)*nstbc_k(-1);
end;

initval;
  nsfc1=1;
  nsfc2=1;
  nsbc1=1;
  nsbc2=1;
  nsbc3=1;
  nstbs=1;
  nsfs=1;
  nsbs=1;
  nstbc_y = 1.1;
  nstbc_c = 0.8;
  nstbc_h = 0.3;
  nstbc_k = 11.1;
end;

steady(solve_algo = @{solve_algo}, tolf = 1e-11);

@#if block
model_info;
@#endif

@#if bytecode
print_bytecode_static_model;
print_bytecode_dynamic_model;
@#endif

check;

shocks;
  var e_R;
  periods 3;
  values 0.01;

  var e_q;
  periods 1;
  values 0.5;

  var e_eb;
  periods 19;
  values 1;

  var e_sfc1;
  periods 9;
  values 0.2;

  var e_sfc2;
  periods 12;
  values -0.4;

  var e_sbc;
  periods 17 18;
  values 0.3 -0.1;

  var e_stbs;
  periods 15;
  values -0.15;

  var e_sbs;
  periods 5;
  values 0.2;

  var e_sfs;
  periods 7;
  values -0.2;

  var e_nsfc1;
  periods 9;
  values 0.2;

  var e_nsfc2;
  periods 12;
  values -0.4;

  var e_nsbc;
  periods 17 18;
  values 0.3 -0.1;

  var e_nstbs;
  periods 15;
  values -0.15;

  var e_nsbs;
  periods 5;
  values 0.2;

  var e_nsfs;
  periods 7;
  values -0.2;

  var e_nstbc;
  periods 3;
  values 0.009;
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver(no_homotopy, markowitz=0, stack_solve_algo = @{stack_solve_algo});
