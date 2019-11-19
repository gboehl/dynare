var dx dy;
varexo e_x e_y;

parameters rho_x rho_y;

rho_x = 0.5;
rho_y = -0.3;

model;
dx = rho_x*dx(-1)+e_x;
dy = rho_y*dy(-1)+e_y;
end;

estimated_params;
rho_x,NORMAL_PDF,0.5,0.1;
rho_y,NORMAL_PDF,-0.3,0.1;
stderr e_x,INV_GAMMA_PDF,0.01,inf;
stderr e_y,INV_GAMMA_PDF,0.01,inf;
end;

varobs dx dy;
check;
if isoctave || ~matlab_ver_less_than('7.14') % xlsread is able to read XLSX without Excel installed since R2012a
estimation(datafile='test.xlsx',nobs=1000,mh_replic=2000,mh_jscale=1.3);
else
estimation(datafile='test.xls',nobs=1000,mh_replic=2000,mh_jscale=1.3);
end
