function plot_regimes(regimes,M_,options_)

nperiods = size(regimes,2);
nconstr = length(fieldnames(regimes(1)))/2;
if nconstr ==1
    regime(1) = {'regime'};
    regimestart(1) = {'regimestart'};
else
    regime(1) = {'regime1'};
    regimestart(1) = {'regimestart1'};
    regime(2) = {'regime2'};
    regimestart(2) = {'regimestart2'};
end    

GraphDirectoryName = CheckPath('graphs',M_.dname);

fhandle = dyn_figure(options_.nodisplay,'Name',[M_.fname ' occbin regimes']);


for k=1:nconstr
    subplot(nconstr,1,k)
    for t=1:nperiods
        nregimes_in_t = length(regimes(t).(regime{k}));
        start_periods = regimes(t).(regimestart{k});
        start_periods = [start_periods max(start_periods)];
        for r=1:nregimes_in_t
            isconstrained = regimes(t).(regime{k})(r);
            if isconstrained
                plot(t,start_periods(r),'*r')
                hold all,
                plot([t t],start_periods(r:r+1),'-r')
            else
                plot(t,start_periods(r),'ob')
                hold all,
                plot([t t],start_periods(r:r+1),'-b')
            end
        end
    end
    title(['regime ' int2str(k)])
    xlabel('historic period')
    ylabel('regime expected start')
end
annotation('textbox',[.25,0,.15,.05],'String','Unbinding','Color','blue');
annotation('textbox',[.65,0,.15,.05],'String','Binding','Color','red');


dyn_saveas(fhandle,[GraphDirectoryName, filesep, M_.fname '_occbin_regimes'],options_.nodisplay,options_.graph_format);
