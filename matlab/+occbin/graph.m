function graph(M_, options_, options_occbin_, oo_, var_list)
% function graph(M_, options_, options_occbin_, oo_, var_list)
%
% Inputs:
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - options_occbin_     [structure]     Matlab's structure containing Occbin options
% - oo_                 [structure]     Matlab's structure containing the results
% - var_list            [char]          list of the variables to plot 

% Copyright © 2021-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

options_ = occbin.set_option(options_,options_occbin_,'graph.steady_state');

if ~exist([M_.dname '/graphs'],'dir')
    mkdir(M_.dname,'graphs');
end

fidTeX=write_latex_header(M_,options_);

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

%get endogenous variables
[i_var, number_of_plots_to_draw_endo, index_uniques] = varlist_indices(var_list, M_.endo_names, 1);
var_list_plots=var_list(index_uniques);
var_list_TeX = M_.endo_names_tex(i_var);

data_to_plot(:,:,1)=oo_.occbin.simul.piecewise(:,i_var);
if isfield(oo_.occbin,'simul') && isfield(oo_.occbin.simul,'linear')
    data_to_plot(:,:,2)=oo_.occbin.simul.linear(:,i_var);
    legend_list = {'Piecewise Linear','Linear'};
else
    legend_list = {'Piecewise Linear'};
end

nperiods=size(data_to_plot,1);
ndim=size(data_to_plot,3);

if ~options_.occbin.graph.steady_state
    data_to_plot=data_to_plot-repmat(oo_.occbin.simul.ys(i_var)',nperiods,1,ndim);
end

%get exogenous variables
[i_var_exo, number_of_plots_to_draw_exo, index_uniques] = varlist_indices(var_list, M_.exo_names, 1);
var_list_plots=[var_list_plots; var_list(index_uniques)];
var_list_TeX = [var_list_TeX; M_.exo_names_tex(i_var_exo)];

if number_of_plots_to_draw_exo>0
    exo_index=NaN(number_of_plots_to_draw_exo);
    for ii=1:length(i_var_exo)
        temp_index=find(oo_.occbin.simul.exo_pos==i_var_exo(ii));
        if ~isempty(temp_index)
            exo_index(ii)=temp_index;
        else
            error('%s was not part of the shocks for Occbin.', var_list{i_var_exo(ii)});
        end
    end
    data_to_plot(:,end+1:end+number_of_plots_to_draw_exo,1)=[oo_.occbin.simul.shocks_sequence(:,exo_index); zeros(nperiods-size(oo_.occbin.simul.shocks_sequence,1),number_of_plots_to_draw_exo)];
    data_to_plot(:,end+1:end+number_of_plots_to_draw_exo,2)=NaN;
end

[nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw_endo+number_of_plots_to_draw_exo);

for fig = 1:nbplt
    hh_fig = dyn_figure(options_.nodisplay,'Name',['Occbin simulated paths, figure ' int2str(fig)]);
    for plt = 1:nstar
        if fig==nbplt && ~lr==0
            subplot(lr,lc,plt);
        else
            subplot(nr,nc,plt);
        end
        h_zero=plot([1 nperiods],[0 0],'--k','linewidth',0.5);
        hold on
        h1=plot(1:nperiods,data_to_plot(:,(fig-1)*nstar+plt,1),'b-','linewidth',2); 
        if ndim==2 && (fig-1)*nstar+plt<=number_of_plots_to_draw_endo
            h2=plot(1:nperiods,data_to_plot(:,(fig-1)*nstar+plt,2),'r--','linewidth',2); hold on
        end        
        hold off
        
        max_y = max(max(data_to_plot(:,(fig-1)*nstar+plt,:)));
        min_y = min(min(data_to_plot(:,(fig-1)*nstar+plt,:)));
        
        y_bottom = min_y - .01*abs(min_y);
        
        y_top = max_y + 0.01*abs(max_y);
        if y_bottom==y_top
            y_top=y_bottom+1;
        end
        axis([1 nperiods y_bottom y_top])
        remove_fractional_xticks
        if plt==1
            if ndim==2
                legend([h1,h2],legend_list,'box','off')
            else
                legend(h1,legend_list,'box','off')
            end        
        end
        if options_.TeX
            title(['$' var_list_TeX{(fig-1)*nstar+plt,:} '$'],'Interpreter','latex');
        else
            title(deblank(var_list_plots((fig-1)*nstar+plt,:)),'Interpreter','none');
        end
        if number_of_plots_to_draw_endo+number_of_plots_to_draw_exo==(fig-1)*nstar+plt
            break
        end
    end
    dyn_saveas(hh_fig,[M_.dname, '/graphs/'  M_.fname '_occbin_' int2str(fig)],options_.nodisplay,options_.graph_format);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        if fig==nbplt
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_occbin_%s}\n',options_.figures.textwidth*min(plt/nc,1),[M_.dname, '/graphs/' M_.fname],int2str(fig));
        else
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_occbin_%s}\n',options_.figures.textwidth*min(plt/lc,1),[M_.dname, '/graphs/' M_.fname],int2str(fig));
        end
        fprintf(fidTeX,'\\caption{Simulated time paths over time. Blue solid line: taking occasionally binding constraint into account. Red dashed line: ignoring constraint.}\n');
        fprintf(fidTeX,'\\label{Fig:occbin:%s}\n', int2str(fig));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'%% End Of TeX file. \n');
    fclose(fidTeX);
end


function fidTeX=write_latex_header(M_,options_)

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([M_.dname, '/graphs/' M_.fname '_occbin.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by occbin.make_chart.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
else
    fidTeX =[];
end
