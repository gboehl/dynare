var y 
    y_s 
    R 
    pie $\pi$ 
    dq  
    pie_s 
    de 
    A   
    y_obs   ${y^{obs}}$ 
    pie_obs ${\pi^{obs}}$
    R_obs ${R^{obs}}$;

varexo e_R ${\varepsilon^R}$ 
    e_q ${\varepsilon^q}$ 
    e_ys ${\varepsilon^{ys}}$
    e_pies ${\varepsilon^\pi}$
    e_A ${\varepsilon^A}$;

parameters psi1 ${\psi_1}$ 
    psi2 ${\psi_2}$ 
    psi3 ${\psi_3}$
    rho_R ${\rho_R}$  
    tau    ${\tau}$
    alpha ${\alpha}$
    rr 
    k 
    rho_q ${\rho_q}$
    rho_A ${\rho_A}$
    rho_ys ${\rho_{ys}}$
    rho_pies    ${\rho_\pi}$;

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


model(linear);
y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+k*alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
pie = de+(1-alpha)*dq+pie_s;
R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;
dq = rho_q*dq(-1)+e_q;
y_s = rho_ys*y_s(-1)+e_ys;
pie_s = rho_pies*pie_s(-1)+e_pies;
A = rho_A*A(-1)+e_A;
y_obs = y-y(-1)+A;
pie_obs = 4*pie;
R_obs = 4*R;
end;

shocks;
	var e_R = 1.25^2;
	var e_q = 2.5^2;
	var e_A = 1.89;
	var e_ys = 1.89;
	var e_pies = 1.89;
end;

varobs y_obs R_obs pie_obs dq de;

options_.TeX=1;