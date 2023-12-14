function  plot_irfs(M_,oo_,options_,irf3,irf4)

shocknames = options_.occbin.plot_irf.exo_names;
simulname = options_.occbin.plot_irf.simulname;
if isempty(simulname)
   simulname_ = simulname;
else 
   simulname_ = [ simulname '_' ];
end
vars_irf = options_.occbin.plot_irf.endo_names;
endo_names_long = options_.occbin.plot_irf.endo_names_long;
endo_scaling_factor = options_.occbin.plot_irf.endo_scaling_factor;
length_irf = options_.occbin.plot_irf.tplot;
if isempty(length_irf)
    length_irf = options_.irf;
end

irflocation_lin   = oo_.occbin.linear_irfs;
irflocation_piece = oo_.occbin.irfs;


steps_irf = 1; 
warning('off','all')

DirectoryName = CheckPath('graphs',M_.dname);

iexo=[];
for i=1:size(shocknames,1)
    itemp = strmatch(shocknames{i},M_.exo_names,'exact');
    if isempty(itemp)
        error(['Shock ',shocknames{i},' is not defined!'])
    else
        iexo=[iexo, itemp];
    end
end

ncols     = options_.occbin.plot_irf.ncols;
nrows     = options_.occbin.plot_irf.nrows;
npan      = ncols*nrows;

plot_grid = options_.occbin.plot_irf.grid;
shocksigns = options_.occbin.plot_irf.shocksigns; 
threshold = options_.occbin.plot_irf.threshold; 

% Add steady_state
if options_.occbin.plot_irf.add_steadystate
    add_stst = options_.occbin.plot_irf.add_steadystate;
else
    add_stst = 0;
end 
for sss = 1:numel(shocksigns)
    
    shocksign = shocksigns{sss};
    
    for j=1:size(shocknames,1)
        %shocknames = M_.exo_names{j};

        j1   = 0;
        isub = 0;
        ifig = 0;

        % Variables
        % ----------------------
        for i = 1:length(vars_irf)

            j1=j1+1;
            if mod(j1,npan)==1
               % vector corresponds to [left bottom width height]. 680 and 678 for the left and bottom elements correspond to the default values used by MATLAB while creating a figure and width, .
                hfig = dyn_figure(options_.nodisplay,'name',['OccbinIRFs ' shocknames{j} ' ' simulname ' ' shocksign],'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','portrait','renderermode','auto','position',[10 10 950 650]);
                ifig=ifig+1;
                isub=0;
            end
            isub=isub+1;      
            
            if isempty(endo_scaling_factor)
                exofactor = 1;
            else               
                exofactor = endo_scaling_factor{i};
            end 

            subplot(nrows,ncols,isub)            
            irf_field   = strcat(vars_irf{i,1},'_',shocknames{j},'_',shocksign);
            
            irfvalues   = irflocation_lin.(irf_field);
            if add_stst
                irfvalues = irfvalues + get_mean(vars_irf{i,1});
            end          
            
            irfvalues(abs(irfvalues) <threshold) = 0;
            
            plot(irfvalues(1:steps_irf:length_irf)*exofactor,'linewidth',2);
            hold on

            irfvalues   = irflocation_piece.(irf_field);  
            if add_stst
                irfvalues = irfvalues + get_mean(vars_irf{i,1});
            end
            irfvalues(abs(irfvalues) <threshold) = 0;
            plot(irfvalues(1:steps_irf:length_irf)*exofactor,'r--','linewidth',2);
            
            hold on
            plot(irfvalues(1:steps_irf:length_irf)*0,'k-','linewidth',1.5);
            % Optional additional IRFs
            if nargin > 10
                  irfvalues   = irf3.(irf_field) ;
                  irfvalues(abs(irfvalues) <threshold) = 0;
                  plot(irfvalues(1:steps_irf:length_irf)*exofactor,'k:','linewidth',2);
            end
            if nargin > 11
                  irfvalues   = irf4.(irf_field)   ;   
                  irfvalues(abs(irfvalues) <threshold) = 0;
                  plot(irfvalues(1:steps_irf:length_irf)*exofactor,'g-.','linewidth',2);
            end
      

            if plot_grid
                grid on
            end

            xlim([1 (length_irf/steps_irf)]);
            
            % title
            if isempty(endo_names_long)
                title(regexprep(vars_irf{i},'_',' '))
            else
                title(endo_names_long{i})
            end

            % Annotation Box + save figure
            % ----------------------
            if mod(j1,npan)==0 || (mod(j1,npan)~=0 && i==length(vars_irf))
                annotation('textbox', [0.1,0,0.35,0.05],'String', 'Linear','Color','Blue','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.55,0,0.35,0.05],'String', 'Piecewise','Color','Red','horizontalalignment','center','interpreter','none');
                dyn_saveas(hfig,[DirectoryName,filesep,M_.fname,'_irf_occbin_',simulname_,shocknames{j},'_',shocksign,'_',int2str(ifig)],options_.nodisplay,options_.graph_format);
            end
        end 
    end
end

warning('on','all')
end