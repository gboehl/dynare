function []=graph_decomp(z,shock_names,endo_names,i_var,initial_date,M_,options_)
% []=graph_decomp(z,shock_names,endo_names,i_var,initial_date,M_,options_)
% Plots the results from the shock_decomposition command
%
% Inputs
%   z               [n_var*(nshock+2)*nperiods]     shock decomposition array, see shock_decomposition.m for details
%   shock_names     [endo_nbr*string length]        shock names from M_.exo_names
%   endo_names      [exo_nbr*string length]         variable names from M_.endo_names
%   i_var           [n_var*1]                       vector indices of requested variables in M_.endo_names and z
%   initial_date    [dseries object]                first period of decomposition to plot
%   M_              [structure]                     Dynare model structure
%   options_        [structure]                     Dynare options structure

% Copyright © 2010-2023 Dynare Team
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

if ~options_.plot_shock_decomp.expand
    GraphDirectoryName = CheckPath('graphs',M_.dname);
end
new_colormap = options_.plot_shock_decomp.colormap;

% number of components equals number of shocks + 1 (initial conditions)
comp_nbr = size(z,2)-1;

SteadyState=[];
fig_mode='';
fig_mode1='';
% fig_name='';
% screen_shocks=0;
opts_decomp = options_.plot_shock_decomp;
if isfield(opts_decomp,'steady_state')
    SteadyState = opts_decomp.steady_state;
end
if ~isempty(opts_decomp.type)
    fig_mode = opts_decomp.type;
    fig_mode1 = ['_' fig_mode];
    fig_mode = [fig_mode '_'];
end

if isfield(opts_decomp,'init_cond_decomp')
    init_cond_decomp = opts_decomp.init_cond_decomp ;
else
    init_cond_decomp = 0;
end

fig_name_long = opts_decomp.fig_name;

use_shock_groups = options_.plot_shock_decomp.use_shock_groups;
screen_shocks = opts_decomp.screen_shocks;
if ~isempty(use_shock_groups) || comp_nbr<=18
    screen_shocks=0;
end
if use_shock_groups
    shock_groups = M_.shock_groups.(use_shock_groups);
    shock_ind = fieldnames(shock_groups);
end
if screen_shocks
    fig_name_long = [fig_name_long ' SCREEN'];
end

fig_name=strrep(fig_name_long, '(given ', '');
fig_name=strrep(fig_name, '(vintage ', '');
fig_name=regexprep(fig_name, ' ', '_');
fig_name=strrep(fig_name, '.', '');
fig_name=strrep(fig_name, '-', '');
fig_name=strrep(fig_name, ')', '');
fig_name=strrep(fig_name, '(', '');

interactive = opts_decomp.interactive;


gend = size(z,3);
if isempty(initial_date)
    x = 0:gend;
    freq = 1;
else
    freq = initial_date.freq;
    initial_period = double(initial_date);
    x = initial_period-1/freq:(1/freq):initial_period+(gend-1)/freq;
end

nvar = length(i_var);

%% write LaTeX-Header
if options_.TeX && any(strcmp('eps',cellstr(options_.plot_shock_decomp.graph_format))) && ~options_.plot_shock_decomp.expand
    fidTeX = fopen([GraphDirectoryName, filesep, M_.fname '_shock_decomp' fig_mode1 fig_name '.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by Dynare''s graph_decomp.m.\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

if init_cond_decomp
    preamble_txt = 'Initial condition decomposition';
    preamble_figname = '_initval_decomposition_';
else
    preamble_figname = '_shock_decomposition_';
    if opts_decomp.vintage && opts_decomp.realtime>1
        preamble_txt = 'Shock decomposition';
    else
        preamble_txt = 'Historical shock decomposition';
    end
end

if ~(screen_shocks && comp_nbr>18)
    screen_shocks=0;
end
comp_nbr0=comp_nbr;
%%plot decomposition
for j=1:nvar
    z1 = squeeze(z(i_var(j),:,:));
    if screen_shocks
        [~, isort] = sort(mean(abs(z1(1:end-2,:)')), 'descend');
        labels = char(char(shock_names(isort(1:16))),'Others', 'Initial values');
        zres = sum(z1(isort(17:end),:),1);
        z1 = [z1(isort(1:16),:); zres; z1(comp_nbr0:end,:)];
        comp_nbr=18;
        func = @(x) colorspace('RGB->Lab',x);
        new_colormap = distinguishable_colors(size(z1,1)-1,'w',func);
        new_colormap(end,:) = [0.7 0.7 0.7];
    else
        labels = char(char(shock_names),'Initial values');
    end
    xmin = x(1);
    xmax = x(end)+1/freq;
    ix = z1(1:comp_nbr,:) > 0;
    ymax = max(sum(z1(1:comp_nbr,:).*ix))*1.1;
    ix = z1(1:comp_nbr,:) < 0;
    ymin = min(sum(z1(1:comp_nbr,:).*ix))*1.1;
    if ymax-ymin < 1e-6
        continue
    end
    fhandle = dyn_figure(options_.plot_shock_decomp.nodisplay,'Name',[preamble_txt fig_name_long strrep(fig_mode1, '_', ' ') ': ' endo_names{i_var(j)} '.'], 'PaperPositionMode', 'auto','PaperOrientation','landscape','renderermode','auto');
    screensize = get( groot, 'Screensize' );
    set(fhandle,'OuterPosition' ,[50 50 min(1500,screensize(3)-50) min(750,screensize(4)-50)])
    ax=axes('Position',[0.1 0.1 0.6 0.8],'box','on');
    %     plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
    %     axis(ax,[xmin xmax ymin ymax]);
    if strcmp('aoa',options_.plot_shock_decomp.type)
        bgap = 0.15;
    else
        bgap = 0;
    end
    hold on;
    for i=1:gend
        i_1 = i-1;
        yp = 0;
        ym = 0;
        for k = 1:comp_nbr
            zz = z1(k,i);
            if zz > 0
                fill([x(i)+bgap x(i)+bgap x(i+1)-bgap x(i+1)-bgap]+(1/freq/2),[yp yp+zz yp+zz yp],k);
                yp = yp+zz;
            else
                fill([x(i)+bgap x(i)+bgap x(i+1)-bgap x(i+1)-bgap]+(1/freq/2),[ym ym+zz ym+zz ym],k);
                ym = ym+zz;
            end
            hold on;
        end
    end
    plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
    if ~isempty(SteadyState)
        plot(ax,[xmin xmax],[0 0],'--','linewidth',1,'color',[0.7 0.7 0.7])
        if ymin+SteadyState(i_var(j))<0 && ymax+SteadyState(i_var(j))>0
            plot(ax,[xmin xmax],SteadyState(i_var(j))*[-1 -1],'k--','linewidth',1)
            ytick=get(ax,'ytick');
            ytick1=ytick-SteadyState(i_var(j));
            ind1=min(find(ytick1>=ymin));
            ind2=max(find(ytick1<=ymax));
            dytick=ytick(2)-ytick(1);
            if ind1>1
                ytick1  = [ytick1(ind1:end) ytick1(end)+dytick:dytick:ymax];
            elseif ind2<length(ytick)
                ytick1= [sort(ytick1(1)-dytick:-dytick:ymin) ytick1(1:ind2)];
            end
            set(ax,'ytick',ytick1),
        else
            ytick1=get(ax,'ytick');
        end
        ylabel = ytick1'+SteadyState(i_var(j));
        ylabel(abs(ylabel)<eps)=0;
        set(ax,'yticklabel',num2str(ylabel,'%4.2g'))
    end
    set(ax,'xlim',[xmin xmax]);
    hold off;

    axes('Position',[0.75 0.1 0.2 0.8]);
    axis([0 1 0 1]);
    axis off;
    hold on;
    y1 = 0;
    height = 1/comp_nbr;

    for i=comp_nbr:-1:1
        %     for i=1:comp_nbr
        hl = fill([0 0 0.2 0.2],[y1 y1+0.7*height y1+0.7*height y1],i);
        hold on
        ht = text(0.3,y1+0.3*height,labels(i,:),'Interpreter','none');
        hold on
        if interactive && (~isoctave && ~isempty(use_shock_groups))
            mydata.fig_name = options_.plot_shock_decomp.fig_name(2:end);
            mydata.use_shock_groups = options_.plot_shock_decomp.use_shock_groups;
            mydata.shock_group = shock_groups.(shock_ind{i});
            mydata.shock_decomp = options_.shock_decomp;
            mydata.plot_shock_decomp = options_.plot_shock_decomp;
            mydata.first_obs = options_.first_obs;
            mydata.nobs = options_.nobs;
            mydata.plot_shock_decomp.zfull = options_.plot_shock_decomp.zfull(i_var(j),:,:);
            mydata.endo_names = endo_names(i_var(j));
            mydata.endo_names_tex = M_.endo_names_tex(i_var(j));
            mydata.exo_names = M_.exo_names;
            if ~isempty(mydata.shock_group.shocks)
                c = uicontextmenu;
                hl.UIContextMenu=c;
                browse_menu = uimenu(c,'Label','Browse group');
                expand_menu = uimenu(c,'Label','Expand group','Callback',['expand_group(''' mydata.plot_shock_decomp.use_shock_groups ''',''' mydata.plot_shock_decomp.orig_varlist{j} ''',' int2str(i) ')']);
                set(expand_menu,'UserData',mydata,'Tag',['group' int2str(i)]);
                save_expand2xls_menu = uimenu(c,'Label','Export group to xls','Callback',['expand_group(''' mydata.plot_shock_decomp.use_shock_groups ''',''' mydata.plot_shock_decomp.orig_varlist{j} ''',' int2str(i) ', 1)']);
                set(save_expand2xls_menu,'Tag',['xls_group' int2str(i)]);
                for jmember = mydata.shock_group.shocks
                    uimenu('parent',browse_menu,'Label',char(jmember))
                end
                ht.UIContextMenu=c;
            end
        end
        y1 = y1 + height;
    end

    if ~isempty(new_colormap)
        colormap(new_colormap)
    end
    hold off
    if ~options_.plot_shock_decomp.expand

        dyn_saveas(fhandle,[GraphDirectoryName, filesep, M_.fname,preamble_figname,endo_names{i_var(j)},fig_mode1,fig_name],options_.plot_shock_decomp.nodisplay,options_.plot_shock_decomp.graph_format);
        if options_.TeX && any(strcmp('eps',cellstr(options_.plot_shock_decomp.graph_format)))
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s/graphs/%s%s}\n',M_.fname,M_.fname,[preamble_figname endo_names{i_var(j)} fig_mode1 fig_name]);
            fprintf(fidTeX,'\\label{Fig:shock_decomp:%s}\n',[fig_mode endo_names{i_var(j)} fig_name]);
            fprintf(fidTeX,['\\caption{' preamble_txt fig_name_long strrep(fig_mode1, '_',  ' ') ': $ %s $.}\n'],M_.endo_names_tex{i_var(j)});
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end
    else
        if ~isempty(options_.plot_shock_decomp.filepath)
            dyn_saveas(fhandle,[options_.plot_shock_decomp.filepath, filesep, M_.fname,preamble_figname,endo_names{i_var(j)},fig_mode1,fig_name],options_.plot_shock_decomp.nodisplay,options_.plot_shock_decomp.graph_format);
        end
    end

end

%% write LaTeX-Footer
if options_.TeX && any(strcmp('eps',cellstr(options_.plot_shock_decomp.graph_format))) &&  ~options_.plot_shock_decomp.expand
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
end