function graph_comparison_irfs(matched_irfs,irf_model_varobs,varobs_id,irf_horizon,relative_irf,endo_names,endo_names_tex,exo_names,exo_names_tex,dname,fname,graph_format,TeX,nodisplay,figures_textwidth)
% graph_comparison_irfs(matched_irfs,irf_model_varobs,varobs_id,irf_horizon,relative_irf,endo_names,endo_names_tex,exo_names,exo_names_tex,dname,fname,graph_format,TeX,nodisplay,figures_textwidth)
% -------------------------------------------------------------------------
% Plots and saves to disk the comparison of the selected data IRFs and corresponding model IRfs
% -------------------------------------------------------------------------
% INPUTS
% matched_irfs:      [matrix]   information on matched data IRFs
% irf_model_varobs:  [matrix]   model IRFs for observable variables
% varobs_id:         [vector]   index for observable variables in endo_names
% irf_horizon:       [scalar]   maximum horizon of IRFs
% relative_irf:      [boolean]  if true, plots normalized IRFs
% endo_names:        [cell]     names of endogenous variables
% endo_names_tex:    [cell]     names of endogenous variables in latex
% exo_names:         [cell]     names of exogenous variables
% exo_names_tex:     [cell]     names of exogenous variables in latex
% dname:             [string]   name of the directory where to save the graphs
% fname:             [string]   name of the mod file
% graph_format:      [cell]     format of the graphs
% TeX:               [boolean]  if true, uses latex for plots
% nodisplay:         [boolean]  if true, does not display the graphs
% figures_textwidth: [scalar]   textwidth used in plots
% -------------------------------------------------------------------------
% OUTPUT
% No output, just displays and saves to disk the graphs
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%  o dyn_figure
%  o dyn_saveas
%  o remove_fractional_xticks
%  o CheckPath
%  o pltorg
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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


graph_directory_name = CheckPath('graphs',dname);
latex_directory_name = CheckPath('latex',dname);
if TeX && any(strcmp('eps',cellstr(graph_format)))
    fid_TeX = fopen([latex_directory_name '/' fname '_irf_matching_plot.tex'],'w');
    fprintf(fid_TeX,'%% TeX eps-loader file generated by mom.run.m (Dynare).\n');
    fprintf(fid_TeX,['%% ' datestr(now,0) '\n']);
    fprintf(fid_TeX,' \n');
end
unique_shock_entries = unique(matched_irfs(:, 2));
colDarkGrey = [0.3, 0.3, 0.3]; % dark grey
for jexo = unique_shock_entries' % loop over cell with shock names
    unique_variables = unique(matched_irfs(ismember(matched_irfs(:, 2),jexo), 1));
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(length(unique_variables));
    fig = 0;
    for jvar = 1:length(unique_variables)
        % get data points, note that periods and values can span over multiple rows 
        jj = ismember(matched_irfs(:,1), unique_variables(jvar)) & ismember(matched_irfs(:,2), jexo);
        IRF_PERIODS = []; IRF_VALUES = [];
        for kk = 1:size(matched_irfs{jj,3},1)
            irf_periods = matched_irfs{jj,3}{kk,1};            
            irf_values = matched_irfs{jj,3}{kk,2};
            if length(irf_values)==1
                irf_values = repmat(irf_values,length(irf_periods),1);
            end            
            IRF_PERIODS = [IRF_PERIODS; irf_periods(:)];
            IRF_VALUES = [IRF_VALUES; irf_values(:)];
        end
        
        if jvar==1 || ~( (fig-1)*nstar<jvar && jvar<=fig*nstar )
            fig = fig+1;
            fig_irf = dyn_figure(nodisplay,'Name',['IRF matching shock to ' jexo{:} ' figure ' int2str(fig)]);
        end
        plt = jvar-(fig-1)*nstar;        
        if nbplt>1 && fig==nbplt
            subplot(lr,lc,plt);
        else
            subplot(nr,nc,plt);
        end
        plt_data = plot(IRF_PERIODS,IRF_VALUES,'h', 'MarkerEdgeColor',colDarkGrey,'MarkerFaceColor',colDarkGrey,'MarkerSize',8);
        hold on
        plt_model = plot(1:irf_horizon, irf_model_varobs(:,varobs_id==find(ismember(endo_names,unique_variables(jvar))) , ismember(exo_names,jexo)),'-k','linewidth',2);
        hold on
        plot([1 irf_horizon],[0 0],'-r','linewidth',1);
        hold off
        xlim([1 irf_horizon]);
        remove_fractional_xticks
        if TeX
            title(['$' endo_names_tex{ismember(endo_names,unique_variables(jvar))} '$'],'Interpreter','latex');
        else
            title(unique_variables{jvar},'Interpreter','none');
        end
        set(gca,'FontSize',12);
        if (plt==nstar) || jvar==length(unique_variables)
            % Adding a legend at the bottom
            axes('Position',[0, 0, 1, 1],'Visible','off');
            lgd = legend([plt_data,plt_model],{'Data', 'Model'}, 'Location', 'southeast','NumColumns',2,'FontSize',14);
            if ~isoctave
                lgd.Position = [0.37 0.01 lgd.Position(3) lgd.Position(4)];
            end
            dyn_saveas(fig_irf,[graph_directory_name filesep fname '_matched_irf_' jexo{:} int2str(fig)],nodisplay,graph_format);
            if TeX && any(strcmp('eps',cellstr(graph_format)))
                fprintf(fid_TeX,'\\begin{figure}[H]\n');
                fprintf(fid_TeX,'\\centering \n');
                fprintf(fid_TeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_matched_irf_%s%s}\n',figures_textwidth*min(plt/nc,1),[graph_directory_name '/' fname],jexo{:},int2str(fig));
                if relative_irf
                    fprintf(fid_TeX,'\\caption{Relative impulse response functions (orthogonalized shock to $%s$).}', jexo{:});
                else
                    fprintf(fid_TeX,'\\caption{Impulse response functions (orthogonalized shock to $%s$).}', jexo{:});
                end
                fprintf(fid_TeX,'\\label{Fig:MatchedIRF:%s:%s}\n', jexo{:},int2str(fig));
                fprintf(fid_TeX,'\\end{figure}\n');
                fprintf(fid_TeX,' \n');
            end
        end
    end
end