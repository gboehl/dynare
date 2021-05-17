function make_chart(titlelist,legendlist,figlabel,ylabels,data_series)
% function make_chart(titlelist,legendlist,figlabel,ylabels,zdata)

% Original authors: Luca Guerrieri and Matteo Iacoviello
% Original file downloaded from:
% https://www.matteoiacoviello.com/research_files/occbin_20140630.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Luca Guerrieri and Matteo Iacoviello (2015): "OccBin: A toolkit for solving
% dynamic models with occasionally binding constraints easily"
% Journal of Monetary Economics 70, 22-38

titlelist = char(strrep(cellstr(titlelist),'_','.'));

ndsets=size(data_series,3);       % default, changed below as applicable
nperiods = size(data_series,1);

xvalues = (1:nperiods)';
nvars = size(titlelist,1);

[nbplt,nr,nc,lr,lc,nstar] = pltorg(nvars);

style_cell={'k','r','b','g','c','y'};

for fig = 1:nbplt
    figure('Name',[figlabel, ', Figure ' int2str(fig)]);
    for plt = 1:nstar
        h1=NaN(ndsets);
        if fig==nbplt && ~lr==0
            subplot(lr,lc,plt);
        else
            subplot(nr,nc,plt);
        end        
        for data_set_iter=1:ndsets
            h1(data_set_iter)=plot(xvalues,data_series(:,(fig-1)*nstar+plt,data_set_iter),style_cell{1+mod(data_set_iter-1,length(style_cell))},'linewidth',2);
            hold on
        end
        grid on
            
        max_y = max(max(data_series(:,(fig-1)*nstar+plt,:)));
        min_y = min(min(data_series(:,(fig-1)*nstar+plt,:)));
        
        y_bottom = min_y - .01*abs(min_y);
        
        y_top = max_y + 0.01*abs(max_y);
        if y_bottom==y_top
            y_top=y_bottom+1;
        end
        
        axis([1 nperiods y_bottom y_top])
        if plt==1
            if numel(strvcat(legendlist(1,:)))
                h=legend(legendlist,'Location','Northwest','Fontsize',8);
            end
        end
        
        title(titlelist(plt,:),'Fontsize',11);
        ylabel(ylabels(plt,:))
        if nvars==(fig-1)*nstar+plt
            break
        end
    end
end