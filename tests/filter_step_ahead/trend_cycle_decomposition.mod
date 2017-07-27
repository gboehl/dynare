var y yp z mu;

varexo ez eyp emu;

parameters alpha;

alpha = .889;

model(linear);
  y  = yp + z;
  yp = mu + yp(-1) + eyp;
  mu = mu(-1) + emu;
  z  = alpha*z(-1) + ez;
end;

initval;
y=8.655680;
z=0;
yp=8.655680;
mu=0;
end;

steady(nocheck);

estimated_params;
stderr emu     , inv_gamma_pdf,  0.002  , inf;
stderr eyp     , inv_gamma_pdf,  0.002  , inf;
stderr ez        , inv_gamma_pdf,   0.06 , inf;
alpha, normal_pdf, 0.9, 0.1;
end;

varobs y; 

estimation(datafile=trend_cycle_decomposition_data,nobs=82, mh_replic=2000, mode_compute=4, mh_nblocks=1, mh_jscale=0.3, filtered_vars, smoother, diffuse_filter) yp z; 
