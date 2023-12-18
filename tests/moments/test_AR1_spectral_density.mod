var white_noise ar1;
varexo e;

parameters phi;

phi=0.9;

model;
white_noise=e;
ar1=phi*ar1(-1)+e;

end;

shocks;
var e = 1;
end;

options_.SpectralDensity.trigger=1;

options_.bandpass.indicator=0;

stoch_simul(order=1,nofunctions,hp_filter=0,irf=0,periods=1000000,filtered_theoretical_moments_grid=2048);
send_endogenous_variables_to_workspace;

white_noise_sample=white_noise;

theoretical_spectrum_white_noise=1^2/(2*pi); %Hamilton (1994), 6.1.9
if max(abs(oo_.SpectralDensity.density(strmatch('white_noise',M_.endo_names,'exact'),:)-theoretical_spectrum_white_noise))>1e-10
   error('Spectral Density is wrong') 
end

theoretical_spectrum_AR1=1/(2*pi)*(1^2./(1+phi^2-2*phi*cos(oo_.SpectralDensity.freqs))); %Hamilton (1994), 6.1.13
if max(abs(oo_.SpectralDensity.density(strmatch('ar1',M_.endo_names,'exact'),:)-theoretical_spectrum_AR1'))>1e-10
   error('Spectral Density is wrong') 
end

stoch_simul(order=1,nofunctions,hp_filter=1600,irf=0,periods=0);
lambda=options_.hp_filter;
Kalman_gain=(4*lambda*(1 - cos(oo_.SpectralDensity.freqs)).^2 ./ (1 + 4*lambda*(1 - cos(oo_.SpectralDensity.freqs)).^2));
theoretical_spectrum_white_noise_hp_filtered=1^2/(2*pi)*Kalman_gain.^2; %Hamilton (1994), 6.1.9
if max(abs(oo_.SpectralDensity.density(strmatch('white_noise',M_.endo_names,'exact'),:)-theoretical_spectrum_white_noise_hp_filtered'))>1e-10
   error('Spectral Density is wrong') 
end

theoretical_spectrum_AR1_hp_filtered=1/(2*pi)*(1^2./(1+phi^2-2*phi*cos(oo_.SpectralDensity.freqs))).*Kalman_gain.^2; %Hamilton (1994), 6.1.13
if max(abs(oo_.SpectralDensity.density(strmatch('ar1',M_.endo_names,'exact'),:)-theoretical_spectrum_AR1_hp_filtered'))>1e-10
   error('Spectral Density is wrong') 
end

options_.hp_filter=0;
stoch_simul(order=1,nofunctions,bandpass_filter=[6 32],irf=0);

theoretical_spectrum_white_noise=repmat(theoretical_spectrum_white_noise,1,options_.filtered_theoretical_moments_grid);
passband=oo_.SpectralDensity.freqs>=2*pi/options_.bandpass.passband(2) & oo_.SpectralDensity.freqs<=2*pi/options_.bandpass.passband(1); 
if max(abs(oo_.SpectralDensity.density(strmatch('white_noise',M_.endo_names,'exact'),passband)-theoretical_spectrum_white_noise(passband)))>1e-10
   error('Spectral Density is wrong') 
end
if max(abs(oo_.SpectralDensity.density(strmatch('white_noise',M_.endo_names,'exact'),~passband)-0))>1e-10
   error('Spectral Density is wrong') 
end

if max(abs(oo_.SpectralDensity.density(strmatch('ar1',M_.endo_names,'exact'),passband)-theoretical_spectrum_AR1(passband)'))>1e-10
   error('Spectral Density is wrong') 
end
if max(abs(oo_.SpectralDensity.density(strmatch('ar1',M_.endo_names,'exact'),~passband)-0))>1e-10
   error('Spectral Density is wrong') 
end


% [pow,f]=psd(a_sample,1024,1,[],512);
% figure
% plot(f,pow/(2*pi))
% 
% % figure
% % [pow,f]=psd(a_sample,1000,1,[],500);
% % plot(f(3:end)*2*pi,pow(3:end)/(2*pi));
