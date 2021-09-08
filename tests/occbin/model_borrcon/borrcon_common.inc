steady;
check;

shocks(surprise,overwrite);
var u;
periods 10, 30;
values 0.03, -0.03;
end;

occbin_setup(simul_periods=80);
occbin_solver(simul_maxit=11,simul_curb_retrench);
% 
titlelist = char('c (consumption)','b (borrowing)','y (income)','lb (multiplier)');
percent = 'Percent';
level = 'Level';
ylabels = char(percent,percent,percent,level);
figtitle = 'Simulated variables';
legendlist = cellstr(char('Piecewise Linear','Linear'));

if isoctave || ~matlab_ver_less_than('8.4')
    % For some unknown reason, MATLAB R2014a crashes with TeX labels
    % when -nodisplay is passed (as in the testsuite)
    options_.TeX=1;
end
occbin_graph;
occbin_graph c_hat b_hat y_hat;

write_latex_original_model(write_equation_tags);
//collect_latex_files;

oo_= occbin.unpack_simulations(M_,oo_,options_);

line1=100*[(oo_.occbin.endo_piecewise.c-oo_.occbin.endo_ss.c)/oo_.occbin.endo_ss.c, ...
        (oo_.occbin.endo_piecewise.b-oo_.occbin.endo_ss.b)/oo_.occbin.endo_ss.b, ...
        (oo_.occbin.endo_piecewise.y-oo_.occbin.endo_ss.y)/oo_.occbin.endo_ss.y, ...
        oo_.occbin.endo_piecewise.lb/100];
line2=100*[(oo_.occbin.endo_linear.c-oo_.occbin.endo_ss.c)/oo_.occbin.endo_ss.c, ...
        (oo_.occbin.endo_linear.b-oo_.occbin.endo_ss.b)/oo_.occbin.endo_ss.b, ...
        (oo_.occbin.endo_linear.y-oo_.occbin.endo_ss.y)/oo_.occbin.endo_ss.y, ...
        oo_.occbin.endo_linear.lb/100];

occbin.make_chart(titlelist,legendlist,figtitle,ylabels,cat(3,line1,line2));