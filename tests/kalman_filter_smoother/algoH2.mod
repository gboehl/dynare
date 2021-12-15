var z dw dx dy dc1 dc2;
varexo e_w e_x e_y e_z;

parameters rho_w rho_x rho_y rho_z a1 a2 a3 b c;

model(linear);
dw = rho_w*dw(-1)+a1*(dc1(-1))+e_w;
dx = rho_x*dx(-1)+a2*(dc1(-1))+e_x;
dy = rho_y*dy(-1)+a3*(dc2(-1))+e_y;
z = rho_z*z(-1)+dw-dx+e_z;
dc1 = dc1(-1)+dx-b*dy-c*dw;
dc2 = dc2(-1)+dx-b*dy;
end;

estimated_params;
rho_w, normal_pdf, 0.5,0.2;
rho_x, normal_pdf, 0.5,0.2;
rho_y, normal_pdf, 0.5,0.2;
rho_z, normal_pdf, 0.8,0.2;

a1, normal_pdf, 0.1,0.2;
a2, normal_pdf,  -0.1,0.2;
a3, normal_pdf,  0.1,0.2;
b , normal_pdf,  1,0.2;
c , normal_pdf,  1,0.2;

stderr e_w, uniform_pdf,,, 0.01, 0.1;
stderr e_x, uniform_pdf,,, 0.01, 0.1;
stderr e_y, uniform_pdf,,, 0.01, 0.1;
stderr e_z, uniform_pdf,,, 0.01, 0.1;

stderr z, normal_pdf, 0, 0.1;
end;

varobs dw dx dy z;
       
//estimation(datafile=data_algo,first_obs=1000,nobs=200,mh_replic=0,kalman_algo=2);
estimation(datafile=data_algo,first_obs=1000,nobs=200,mh_replic=0,mode_compute=0,mode_file='algoH1/Output/algoH1_mode',kalman_algo=2,filtered_vars,smoothed_state_uncertainty);

//checking smoother consistency
X = oo_.SmoothedVariables;
S = [X.z X.dw X.dx X.dy X.dc1 X.dc2];
X = oo_.SmoothedShocks;
E = [X.e_w X.e_x X.e_y X.e_z];
A = oo_.dr.ghx;
B = oo_.dr.ghu;
err = zeros(6,200);
for t=2:200;
    err(:,t) = S(t,:)'-A*S(t-1,:)'-B*E(t,:)';
end;
if max(max(abs(err))) > 1e-10;
   error('Test fails');
end;    

d=load('data_algo');
dat = [d.dw d.dx d.dy d.z];
X = oo_.SmoothedMeasurementErrors;
ME = [X.dw X.dx X.dy X.z];
if max(max(abs(dat(1000:1199,:)-S(:,[2:4 1])-ME))) > 1e-10;
   error('Test fails');
end;    

o1 = load(['algoH1' filesep 'Output' filesep 'algoH1_results.mat']);
obj_endo={'SmoothedVariables'; 'FilteredVariables'; 'UpdatedVariables'};
obj_exo = {'SmoothedShocks';}; 
nobj_endo = size(obj_endo,1);
nobj_exo = size(obj_exo,1);
for i=1:nobj_endo;
    err_endo = zeros(eval(['size(oo_.' obj_endo{i} '.' M_.endo_names{1} ',1);']),M_.endo_nbr);
    for j=1:M_.endo_nbr;
        var1 = eval(['o1.oo_.' obj_endo{i} '.' M_.endo_names{j}]);
        var2 = eval(['oo_.' obj_endo{i} '.' M_.endo_names{j}]);
        err_endo(:,j) = var1-var2;
    end;
    if max(max(abs(err_endo))) > 1e-10;
       error('Test fails');
    end;    
end;


err_exo = zeros(200,M_.exo_nbr,nobj_exo);
for i=1:nobj_exo;
    err_exo = zeros(size(eval(['oo_.' obj_exo{i} '.' M_.exo_names{1}]),1),M_.exo_nbr);
    for j=1:M_.exo_nbr;
        var1 = eval(['o1.oo_.' obj_exo{i} '.' M_.exo_names{j}]);
        var2 = eval(['oo_.' obj_exo{i} '.' M_.exo_names{j}]);
        err_exo(:,j,i) = var1 - var2;
    end;
    if max(max(abs(err_exo))) > 1e-10;
       error('Test fails');
    end;    
end;

if max(max(abs(o1.oo_.SmoothedMeasurementErrors.z - oo_.SmoothedMeasurementErrors.z))) > 1e-10; 
   error('Test fails');
end;    
