function [ ix2, ilogpo2, ModelName, MetropolisFolder, FirstBlock, FirstLine, npar, NumberOfBlocks, nruns, NewFile, MAX_nruns, d, bayestopt_] = ...
    posterior_sampler_initialization(TargetFun, xparam1, vv, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_, dispString)

% Posterior sampler initialization.
%
% INPUTS
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters.
%   o dataset_              data structure (likelihood based) or data moments/irfs (method of moments)
%   o dataset_info          dataset info structure (likelihood based) or info on weighting matrix (method of moments)
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%   o dispString            string to be displayed in the command window
%
% OUTPUTS
%   o ix2                   [double]   (NumberOfBlocks*npar) vector of starting points for different chains
%   o ilogpo2               [double]   (NumberOfBlocks*1) vector of initial posterior values for different chains
%   o ModelName             [string]    name of the mod-file
%   o MetropolisFolder      [string]    path to the Metropolis subfolder
%   o FirstBlock            [scalar]    number of the first MH chain to be run (not equal to 1 in case of recovery)
%   o FirstLine             [double]   (NumberOfBlocks*1) vector of first draw in each chain (not equal to 1 in case of recovery)
%   o npar                  [scalar]    number of parameters estimated
%   o NumberOfBlocks        [scalar]    Number of MCM chains requested
%   o nruns                 [double]   (NumberOfBlocks*1) number of draws in each chain
%   o NewFile               [scalar]    (NumberOfBlocks*1) vector storing the number
%                                       of the first MH-file to created for each chain when saving additional
%                                       draws
%   o MAX_nruns             [scalar]    maximum number of draws in each MH-file on the harddisk
%   o d                     [double]   (p*p) matrix, Cholesky decomposition of the posterior covariance matrix (at the mode).
%   o bayestopt_            [structure] estimation options structure
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2006-2023 Dynare Team
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

%Initialize outputs
ix2 = [];
ilogpo2 = [];
FirstBlock = [];
FirstLine = [];
NewFile = [];

ModelName = M_.fname;
if ~isempty(M_.bvar)
    ModelName = [ModelName '_bvar'];
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
BaseName = [MetropolisFolder filesep ModelName];

NumberOfBlocks = options_.mh_nblck;
nruns = ones(NumberOfBlocks,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
d = chol(vv);

if ~options_.load_mh_file && ~options_.mh_recover
    % Here we start a new Metropolis-Hastings, previous draws are discarded.
    if NumberOfBlocks > 1
        fprintf('%s: Multiple chains mode.\n',dispString);
    else
        fprintf('%s: One Chain mode.\n',dispString);
    end
    % Delete old mh files if any...
    files = dir([BaseName '_mh*_blck*.mat']);
    if length(files)
        delete([BaseName '_mh*_blck*.mat']);
        fprintf('%s: Old mh-files successfully erased!\n',dispString);
    end
    % Delete old Metropolis log file.
    file = dir([ MetropolisFolder '/metropolis.log']);
    if length(file)
        delete([ MetropolisFolder '/metropolis.log']);
        fprintf('%s: Old metropolis.log file successfully erased!\n',dispString)
        fprintf('%s: Creation of a new metropolis.log file.\n',dispString)
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'w');
    fprintf(fidlog,'%% MH log file (Dynare).\n');
    fprintf(fidlog,['%% ' datestr(now,0) '.\n']);
    fprintf(fidlog,' \n\n');
    fprintf(fidlog,'%% Session 1.\n');
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(NumberOfBlocks) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    if isempty(d)
        Prior = dprior(bayestopt_, options_.prior_trunc);
    end
    if options_.mh_initialize_from_previous_mcmc.status
        PreviousFolder0 = options_.mh_initialize_from_previous_mcmc.directory;
        RecordFile0 = options_.mh_initialize_from_previous_mcmc.record;
        PriorFile0 = options_.mh_initialize_from_previous_mcmc.prior;
        if ~isempty(RecordFile0)
            %% check for proper filesep char in user defined paths
            RecordFile0=strrep(RecordFile0,'\',filesep);
            if isempty(dir(RecordFile0))
                fprintf('%s: Wrong value for mh_initialize_from_previous_mcmc_record option.\n',dispString);
                error('%s: Path to record file is not found!',dispString)
            else
                record0=load(RecordFile0);
            end
            record0=record0.record;
            MetropolisFolder0 = fileparts(RecordFile0);
            PreviousFolder0=fileparts(MetropolisFolder0);
            [~, ModelName0]=fileparts(PreviousFolder0);
        else
            %% check for proper filesep char in user defined paths
            PreviousFolder0=strrep(PreviousFolder0,'\',filesep);
            MetropolisFolder0 = [PreviousFolder0 filesep 'metropolis'];
            [~, ModelName0]=fileparts(PreviousFolder0);
            record0=load_last_mh_history_file(MetropolisFolder0, ModelName0);
        end
        if ~isnan(record0.MCMCConcludedSuccessfully) && ~record0.MCMCConcludedSuccessfully
            error('%s: You are trying to load an MCMC that did not finish successfully. Please use ''mh_recover''!',dispString);
        end
%         mh_files = dir([ MetropolisFolder0 filesep ModelName0 '_mh*.mat']);
%         if ~length(mh_files)
%             error('%s: I cannot find any MH file to load here!',dispString)
%         end
        fprintf('%s: Initializing from past Metropolis-Hastings simulations...\n',dispString);
        fprintf('%s: Past MH path %s\n',dispString,MetropolisFolder0);
        fprintf('%s: Past model name %s\n', dispString, ModelName0);
        fprintf(fidlog,'  Loading initial values from previous MH\n');
        fprintf(fidlog,'  Past MH path: %s\n', MetropolisFolder0 );
        fprintf(fidlog,'  Past model name: %s\n', ModelName0);
        fprintf(fidlog,' \n');
        past_number_of_blocks = record0.Nblck;
        if past_number_of_blocks ~= NumberOfBlocks
            fprintf('%s: The specified number of blocks doesn''t match with the previous number of blocks!\n', dispString);
            fprintf('%s: You declared %u blocks, but the previous number of blocks was %u.\n', dispString, NumberOfBlocks, past_number_of_blocks);
            fprintf('%s: I will run the Metropolis-Hastings with %u block.\n', dispString, past_number_of_blocks);
            NumberOfBlocks = past_number_of_blocks;
            options_.mh_nblck = NumberOfBlocks;
        end
        if ~isempty(PriorFile0)
            if isempty(dir(PriorFile0))
                fprintf('%s: Wrong value for mh_initialize_from_previous_mcmc_prior option.', dispString);
                error('%s: Path to prior file is not found!',dispString);
            else
                bayestopt0 = load(PriorFile0);
            end
        else
            bayestopt0 = load([PreviousFolder0 filesep 'prior' filesep 'definition.mat']);
        end
        [~,IA,IB] = intersect(bayestopt_.name,bayestopt0.bayestopt_.name);
        new_estimated_parameters = ~ismember(bayestopt_.name,bayestopt0.bayestopt_.name);
        ix2 = zeros(NumberOfBlocks,npar);
        ilogpo2 = zeros(NumberOfBlocks,1);
        ilogpo2(:,1) = record0.LastLogPost;
        ix2(:,IA) = record0.LastParameters(:,IB);
        for j=1:NumberOfBlocks
            if not(all(ix2(j,:)' >= mh_bounds.lb) && all(ix2(j,:)' <= mh_bounds.ub))
                new_estimated_parameters = logical(new_estimated_parameters + (ix2(j,:)' < mh_bounds.lb));
                new_estimated_parameters = logical(new_estimated_parameters + (ix2(j,:)' > mh_bounds.ub));
            end
        end
    else
        new_estimated_parameters = true(1,npar);
    end
    % Find initial values for the NumberOfBlocks chains...
    if NumberOfBlocks > 1 || options_.mh_initialize_from_previous_mcmc.status% Case 1: multiple chains
        options_=set_dynare_seed_local_options(options_,'default');
        fprintf(fidlog,'  Initial values of the parameters:\n');
        fprintf('%s: Searching for initial values...\n', dispString);
        if ~options_.mh_initialize_from_previous_mcmc.status
            ix2 = zeros(NumberOfBlocks,npar);
            ilogpo2 = zeros(NumberOfBlocks,1);
        end
        for j=1:NumberOfBlocks
            validate    = false;
            init_iter   = 0;
            trial = 1;
            while ~validate && trial <= 10
                if isempty(d)
                    candidate = Prior.draw();
                else
                    if isfield(options_,'mh_init_scale')
                        if trial==1
                            fprintf('\nposterior_sampler_initialization: the mh_init_scale-option is deprecated. You should use the mh_init_scale_factor-option instead.\n')
                        end
                        candidate = rand_multivariate_normal( transpose(xparam1), d * options_.mh_init_scale, npar);
                    else
                        candidate = rand_multivariate_normal( transpose(xparam1), d * options_.mh_init_scale_factor*options_.mh_jscale, npar);
                    end
                end
                if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
                    ix2(j,new_estimated_parameters) = candidate(new_estimated_parameters);
                    ilogpo2(j) = - feval(TargetFun,ix2(j,:)',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
                    if isfinite(ilogpo2(j)) % log density is finite
                        fprintf(fidlog,['    Blck ' int2str(j) ':\n']);
                        for i=1:length(ix2(1,:))
                            fprintf(fidlog,['      params:' int2str(i) ': ' num2str(ix2(j,i)) '\n']);
                        end
                        fprintf(fidlog,['      logpo2: ' num2str(ilogpo2(j)) '\n']);
                        validate = true;
                    end
                end
                init_iter = init_iter + 1;
                if  ~validate && init_iter>100
                    fprintf('%s: I couldn''t get a valid initial value in 100 trials.\n', dispString);
                    if options_.nointeractive
                        fprintf('%s: I reduce ''mh_init_scale_factor'' by 10 percent:\n', dispString);
                        if isfield(options_,'mh_init_scale')
                           options_.mh_init_scale = .9*options_.mh_init_scale;
                           fprintf('%s: Parameter ''mh_init_scale'' is now equal to %f.\n',dispString,options_.mh_init_scale);
                        else
                            options_.mh_init_scale_factor = .9*options_.mh_init_scale_factor;
                            fprintf('%s: Parameter ''mh_init_scale_factor'' is now equal to %f.\n', dispString,options_.mh_init_scale_factor);
                        end
                        fprintf('%s: Parameter ''mh_init_scale_factor'' is now equal to %f.\n', dispString,options_.mh_init_scale_factor);
                    else
                        if isfield(options_,'mh_init_scale')
                            fprintf('%s: You should reduce mh_init_scale...\n',dispString);
                            fprintf('%s: Parameter ''mh_init_scale'' is equal to %f.\n',dispString,options_.mh_init_scale);
                            options_.mh_init_scale_factor = input(sprintf('%s: Enter a new value...  ',dispString));
                        else
                            fprintf('%s: You should reduce ''mh_init_scale_factor''...\n',dispString);
                            fprintf('%s: Parameter ''mh_init_scale_factor'' is equal to %f.\n',dispString,options_.mh_init_scale_factor);
                            options_.mh_init_scale_factor = input(sprintf('%s: Enter a new value...  ',dispString));
                        end
                    end
                    trial = trial+1;
                end
            end
            if ~validate && trial>10
                fprintf('%s: I''m unable to find a starting value for block %u.', dispString, j);
                fclose(fidlog);
                return
            end
        end
        fprintf(fidlog,' \n');
        fprintf('%s: Initial values found!\n\n',dispString);
    else% Case 2: one chain (we start from the posterior mode)
        fprintf(fidlog,'  Initial values of the parameters:\n');
        candidate = transpose(xparam1(:));%
        if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
            ix2 = candidate;
            ilogpo2 = - feval(TargetFun,ix2',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
            fprintf('%s: Initialization at the posterior mode.\n\n',dispString);
            fprintf(fidlog,['    Blck ' int2str(1) 'params:\n']);
            for i=1:length(ix2(1,:))
                fprintf(fidlog,['      ' int2str(i)  ':' num2str(ix2(1,i)) '\n']);
            end
            fprintf(fidlog,['    Blck ' int2str(1) 'logpo2:' num2str(ilogpo2) '\n']);
        else
            fprintf('%s: Initialization failed...\n',dispString);
            fprintf('%s: The posterior mode lies outside the prior bounds.\n',dispString);
            fclose(fidlog);
            return
        end
        fprintf(fidlog,' \n');
    end
    fprintf(fidlog,' \n');
    FirstBlock = 1;
    FirstLine = ones(NumberOfBlocks,1);
    NewFile = ones(NumberOfBlocks,1);
    % Delete the mh-history files
    delete_mh_history_files(MetropolisFolder, ModelName);
    %  Create a new record structure
    fprintf('%s: Write details about the MCMC... ', dispString);
    AnticipatedNumberOfFiles = ceil(nruns(1)/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = nruns(1) - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.Sampler = options_.posterior_sampler_options.posterior_sampling_method;
    record.Nblck = NumberOfBlocks;
    record.MhDraws = zeros(1,3);
    record.MhDraws(1,1) = nruns(1);
    record.MhDraws(1,2) = AnticipatedNumberOfFiles;
    record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.MAX_nruns=MAX_nruns;
    record.AcceptanceRatio = zeros(1,NumberOfBlocks);
    record.FunctionEvalPerIteration = NaN(1,NumberOfBlocks);
    if options_.mh_initialize_from_previous_mcmc.status
        record.InitialSeeds = record0.LastSeeds;
    else
        for j=1:NumberOfBlocks
            % we set a different seed for the random generator for each block then we record the corresponding random generator state (vector)
            options_=set_dynare_seed_local_options(options_,options_.DynareRandomStreams.seed+j);
            % record.Seeds keeps a vector of the random generator state and not the scalar seed despite its name
            [record.InitialSeeds(j).Unifor,record.InitialSeeds(j).Normal] = get_dynare_random_generator_state();
        end
    end
    record.InitialParameters = ix2;
    record.InitialLogPost = ilogpo2;
    record.LastParameters = zeros(NumberOfBlocks,npar);
    record.LastLogPost = zeros(NumberOfBlocks,1);
    record.LastFileNumber = AnticipatedNumberOfFiles ;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MCMCConcludedSuccessfully = 0;
    record.ProposalScaleVec=bayestopt_.jscale;
    record.ProposalCovariance=d;
    fprintf('Ok!\n');
    id = write_mh_history_file(MetropolisFolder, ModelName, record);
    fprintf('%s: Details about the MCMC are available in %s_mh_history_%u.mat\n\n', dispString,BaseName,id);
    fprintf(fidlog,'  CREATION OF THE MH HISTORY FILE!\n\n');
    fprintf(fidlog,['    Expected number of files per block.......: ' int2str(AnticipatedNumberOfFiles) '.\n']);
    fprintf(fidlog,['    Expected number of lines in the last file: ' int2str(AnticipatedNumberOfLinesInTheLastFile) '.\n']);
    fprintf(fidlog,'\n');
    for j = 1:NumberOfBlocks
        fprintf(fidlog,['    Initial state of the Gaussian random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Normal)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Normal(i)') '\n']);
        end
        fprintf(fidlog,['    Initial state of the Uniform random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Unifor)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Unifor(i)') '\n']);
        end
    end
    fprintf(fidlog,' \n');
    fclose(fidlog);
elseif options_.load_mh_file && ~options_.mh_recover
    % Here we consider previous mh files (previous mh did not crash).
    fprintf('%s: I am loading past Metropolis-Hastings simulations...\n',dispString);
    record=load_last_mh_history_file(MetropolisFolder, ModelName);
    if ~isnan(record.MCMCConcludedSuccessfully) && ~record.MCMCConcludedSuccessfully
        error('%s: You are trying to load an MCMC that did not finish successfully. Please use mh_recover.',dispString);
    end
    record.MCMCConcludedSuccessfully=0; %reset indicator for this run
    mh_files = dir([ MetropolisFolder filesep ModelName '_mh*.mat']);
    if ~length(mh_files)
        error('%s: I cannot find any MH file to load here!',dispString);
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
    fprintf(fidlog,'\n');
    fprintf(fidlog,['%% Session ' int2str(length(record.MhDraws(:,1))+1) '.\n']);
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(NumberOfBlocks) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    past_number_of_blocks = record.Nblck;
    if past_number_of_blocks ~= NumberOfBlocks
        fprintf('%s: The specified number of blocks doesn''t match with the previous number of blocks!\n',dispString);
        fprintf('%s: You declared %u blocks, but the previous number of blocks was %u.\n', dispString,NumberOfBlocks,past_number_of_blocks);
        fprintf('%s: I will run the Metropolis-Hastings with %u blocks.\n', dispString,past_number_of_blocks);
        NumberOfBlocks = past_number_of_blocks;
        options_.mh_nblck = NumberOfBlocks;
    end
    % I read the last line of the last mh-file for initialization of the new Metropolis-Hastings simulations:
    LastFileNumber = record.LastFileNumber;
    LastLineNumber = record.LastLineNumber;
    if LastLineNumber < MAX_nruns
        NewFile = ones(NumberOfBlocks,1)*LastFileNumber;
        FirstLine = ones(NumberOfBlocks,1)*(LastLineNumber+1);
    else
        NewFile = ones(NumberOfBlocks,1)*(LastFileNumber+1);
        FirstLine = ones(NumberOfBlocks,1);
    end
    ilogpo2 = record.LastLogPost;
    ix2 = record.LastParameters;
    [d,bayestopt_,record]=set_proposal_density_to_previous_value(record,options_,bayestopt_,d,dispString);
    FirstBlock = 1;
    NumberOfPreviousSimulations = sum(record.MhDraws(:,1),1);
    fprintf('%s: I am writing a new mh-history file... ',dispString);
    record.MhDraws = [record.MhDraws;zeros(1,3)];
    NumberOfDrawsWrittenInThePastLastFile = MAX_nruns - LastLineNumber;
    NumberOfDrawsToBeSaved = nruns(1) - NumberOfDrawsWrittenInThePastLastFile;
    AnticipatedNumberOfFiles = ceil(NumberOfDrawsToBeSaved/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = NumberOfDrawsToBeSaved - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.LastFileNumber = LastFileNumber + AnticipatedNumberOfFiles;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MhDraws(end,1) = nruns(1);
    record.MhDraws(end,2) = AnticipatedNumberOfFiles;
    record.MhDraws(end,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.MAX_nruns = [record.MAX_nruns;MAX_nruns];
    record.InitialSeeds = record.LastSeeds;
    write_mh_history_file(MetropolisFolder, ModelName, record);
    fprintf('Done.\n')
    if ~options_.use_mh_covariance_matrix
        fprintf('%s: Ok. I have loaded %u simulations.\n\n', dispString,NumberOfPreviousSimulations);
    end
    fclose(fidlog);
elseif options_.mh_recover
    % The previous metropolis-hastings crashed before the end! I try to recover the saved draws...
    fprintf('%s: Recover mode!\n',dispString);
    record=load_last_mh_history_file(MetropolisFolder, ModelName);
    NumberOfBlocks = record.Nblck;% Number of "parallel" mcmc chains.
    options_.mh_nblck = NumberOfBlocks;

    %% check consistency of options
    if record.MhDraws(end,1)~=options_.mh_replic
        fprintf('\n%s: You cannot specify a different mh_replic than in the chain you are trying to recover\n',dispString);
        fprintf('%s: I am resetting mh_replic to %u\n',dispString,record.MhDraws(end,1));
        options_.mh_replic = record.MhDraws(end,1);
        nruns = ones(NumberOfBlocks,1)*options_.mh_replic;
    end

    if ~isnan(record.MAX_nruns(end,1)) %field exists
        if record.MAX_nruns(end,1)~=MAX_nruns
            fprintf('\n%s: You cannot specify a different MaxNumberOfBytes than in the chain you are trying to recover\n',dispString);
            fprintf('%s: I am resetting MAX_nruns to %u\n',dispString,record.MAX_nruns(end,1));
            MAX_nruns=record.MAX_nruns(end,1);
        end
    end
    %% do tentative initialization as if full last run of MCMC were to be redone
    if size(record.MhDraws,1) == 1
        OldMhExists = 0;% The crashed metropolis was the first session.
    else
        OldMhExists = 1;% The crashed metropolis wasn't the first session.
    end
    % Default initialization:
    if OldMhExists
        ilogpo2 = record.LastLogPost;
        ix2 = record.LastParameters;
    else
        ilogpo2 = record.InitialLogPost;
        ix2 = record.InitialParameters;
    end
    % Set NewFile, a NumberOfBlocks*1 vector of integers, and FirstLine (first line), a NumberOfBlocks*1 vector of integers.
    % Relevant for starting at the correct file and potentially filling a file from the previous run that was not yet full
    if OldMhExists
        LastLineNumberInThePreviousMh = record.MhDraws(end-1,3);% Number of lines in the last mh files of the previous session.
        LastFileNumberInThePreviousMh = sum(record.MhDraws(1:end-1,2),1);% Number of mh files in the the previous sessions.
                                                                         %Test if the last mh files of the previous session were not full yet
        if LastLineNumberInThePreviousMh < MAX_nruns%not full
                                                    %store starting point if whole chain needs to be redone
            NewFile = ones(NumberOfBlocks,1)*LastFileNumberInThePreviousMh;
            FirstLine = ones(NumberOfBlocks,1)*(LastLineNumberInThePreviousMh+1);
            LastFileFullIndicator=0;
        else% The last mh files of the previous session were full
            %store starting point if whole chain needs to be redone
            NewFile = ones(NumberOfBlocks,1)*(LastFileNumberInThePreviousMh+1);
            FirstLine = ones(NumberOfBlocks,1);
            LastFileFullIndicator=1;
        end
    else
        LastLineNumberInThePreviousMh = 0;
        LastFileNumberInThePreviousMh = 0;
        NewFile = ones(NumberOfBlocks,1);
        FirstLine = ones(NumberOfBlocks,1);
        LastFileFullIndicator=1;
    end
    if ~isequal(options_.posterior_sampler_options.posterior_sampling_method,'slice')
        [d,bayestopt_,record]=set_proposal_density_to_previous_value(record,options_,bayestopt_,d,dispString);
    end
    %% Now find out what exactly needs to be redone
    % 1. Check if really something needs to be done
    % How many mh files should we have ?
    ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
    ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*NumberOfBlocks;
    % How many mh files do we actually have ?
    AllMhFiles = dir([BaseName '_mh*_blck*.mat']);
    TotalNumberOfMhFiles = length(AllMhFiles)-length(dir([BaseName '_mh_tmp*_blck*.mat']));
    % Quit if no crashed mcmc chain can be found as there are as many files as expected
    if (ExpectedNumberOfMhFilesPerBlock>LastFileNumberInThePreviousMh) && (TotalNumberOfMhFiles==ExpectedNumberOfMhFiles)
        if isnumeric(options_.parallel)
            fprintf('%s: It appears that you don''t need to use the mh_recover option!\n',dispString);
            fprintf('                  You have to edit the mod file and remove the mh_recover option\n');
            fprintf('                  in the estimation command.\n');
            error('%s: mh_recover option not required!',dispString)
        end
    end
    % 2. Something needs to be done; find out what
    % Count the number of saved mh files per block.
    NumberOfMhFilesPerBlock = zeros(NumberOfBlocks,1);
    is_chain_complete = true(NumberOfBlocks,1);
    for b = 1:NumberOfBlocks
        BlckMhFiles = dir([BaseName '_mh*_blck' int2str(b) '.mat']);
        NumberOfMhFilesPerBlock(b) = length(BlckMhFiles)-length(dir([BaseName '_mh_tmp*_blck' int2str(b) '.mat']));
        if (ExpectedNumberOfMhFilesPerBlock==LastFileNumberInThePreviousMh)
            % here we need to check for the number of lines
            tmpdata = load([BaseName '_mh' int2str(LastFileNumberInThePreviousMh) '_blck' int2str(b) '.mat'],'logpo2');
            if length(tmpdata.logpo2) == LastLineNumberInThePreviousMh
                is_chain_complete(b) = false;
            end
        end
    end
    % Find FirstBlock (First block), an integer targeting the crashed mcmc chain.
    FirstBlock = 1; %initialize
    FBlock = zeros(NumberOfBlocks,1);
    while FirstBlock <= NumberOfBlocks
        if  (NumberOfMhFilesPerBlock(FirstBlock) < ExpectedNumberOfMhFilesPerBlock) || not(is_chain_complete(FirstBlock))
            fprintf('%s: Chain %u is not complete!\n', dispString,FirstBlock);
            FBlock(FirstBlock)=1;
        else
            fprintf('%s: Chain %u is complete!\n', dispString,FirstBlock);
        end
        FirstBlock = FirstBlock+1;
    end

    %% 3. Overwrite default settings for
    % How many mh-files are saved in this block?
    ExistingDrawsInLastMCFile=zeros(NumberOfBlocks,1); %initialize: no MCMC draws of current MCMC are in file from last run
    % Check whether last present file is a file included in the last MCMC run

    update_record=0;
    for k=1:NumberOfBlocks
        FirstBlock = k;
        if FBlock(k)
            NumberOfSavedMhFilesInTheCrashedBlck=NumberOfMhFilesPerBlock(k);
            if ~LastFileFullIndicator
                if NumberOfSavedMhFilesInTheCrashedBlck==NewFile(FirstBlock) %only that last file exists, but no files from current MCMC
                    loaded_results=load([BaseName '_mh' int2str(NewFile(FirstBlock)) '_blck' int2str(FirstBlock) '.mat']);
                    %check whether that last file was filled
                    if size(loaded_results.x2,1)==MAX_nruns %file is full
                        NewFile(FirstBlock)=NewFile(FirstBlock)+1; %set first file to be created to next one
                        FirstLine(FirstBlock) = 1; %use first line of next file
                        ExistingDrawsInLastMCFile(FirstBlock)=MAX_nruns-record.MhDraws(end-1,3);
                    else
                        ExistingDrawsInLastMCFile(FirstBlock)=0;
                    end
                end
            elseif LastFileFullIndicator
                ExistingDrawsInLastMCFile(FirstBlock)=0;
                if NumberOfSavedMhFilesInTheCrashedBlck==NewFile(FirstBlock) %only the last file exists, but no files from current MCMC
                    NewFile(FirstBlock)=NewFile(FirstBlock)+1; %set first file to be created to next one
                end
            end
            %     % Correct the number of saved mh files if the crashed Metropolis was not the first session (so
            %     % that NumberOfSavedMhFilesInTheCrashedBlck is the number of saved mh files in the crashed chain
            %     % of the current session).
            %     if OldMhExists
            %         NumberOfSavedMhFilesInTheCrashedBlck = NumberOfSavedMhFilesInTheCrashedBlck - LastFileNumberInThePreviousMh;
            %     end
            %     NumberOfSavedMhFiles = NumberOfSavedMhFilesInTheCrashedBlck+LastFileNumberInThePreviousMh;

            % Correct initial conditions.
            if NumberOfSavedMhFilesInTheCrashedBlck>0 && NumberOfSavedMhFilesInTheCrashedBlck<ExpectedNumberOfMhFilesPerBlock
                loaded_results=load([BaseName '_mh' int2str(NumberOfSavedMhFilesInTheCrashedBlck) '_blck' int2str(FirstBlock) '.mat']);
                ilogpo2(FirstBlock) = loaded_results.logpo2(end);
                ix2(FirstBlock,:) = loaded_results.x2(end,:);
                nruns(FirstBlock)=nruns(FirstBlock)-ExistingDrawsInLastMCFile(FirstBlock)-(NumberOfSavedMhFilesInTheCrashedBlck-LastFileNumberInThePreviousMh)*MAX_nruns;
                %reset seed if possible
                if isfield(loaded_results,'LastSeeds')
                    record.InitialSeeds(FirstBlock).Unifor=loaded_results.LastSeeds.(['file' int2str(NumberOfSavedMhFilesInTheCrashedBlck)]).Unifor;
                    record.InitialSeeds(FirstBlock).Normal=loaded_results.LastSeeds.(['file' int2str(NumberOfSavedMhFilesInTheCrashedBlck)]).Normal;
                else
                    fprintf('%s: You are trying to recover a chain generated with an older Dynare version.\n',dispString);
                    fprintf('%s: I am using the default seeds to continue the chain.\n',dispString);
                end
                if update_record
                    update_last_mh_history_file(MetropolisFolder, ModelName, record);
                else
                    write_mh_history_file(MetropolisFolder, ModelName, record);
                end
                update_record=1;
            end
        else
            loaded_results=load([BaseName '_mh' int2str(ExpectedNumberOfMhFilesPerBlock) '_blck' int2str(FirstBlock) '.mat']);
            ilogpo2(FirstBlock) = loaded_results.logpo2(end);
            ix2(FirstBlock,:) = loaded_results.x2(end,:);
            nruns(FirstBlock)=0;
            %reset seed if possible
            if isfield(loaded_results,'LastSeeds')
                record.LastSeeds(FirstBlock).Unifor=loaded_results.LastSeeds.(['file' int2str(ExpectedNumberOfMhFilesPerBlock)]).Unifor;
                record.LastSeeds(FirstBlock).Normal=loaded_results.LastSeeds.(['file' int2str(ExpectedNumberOfMhFilesPerBlock)]).Normal;
            else
                fprintf('%s: You are trying to recover a chain generated with an older Dynare version.\n',dispString);
                fprintf('%s: I am using the default seeds to continue the chain.\n',dispString);
            end
            if isfield(loaded_results,'accepted_draws_this_chain')
                record.AcceptanceRatio(FirstBlock)=loaded_results.accepted_draws_this_chain/record.MhDraws(end,1);
                record.FunctionEvalPerIteration(FirstBlock)=loaded_results.accepted_draws_this_chain/record.MhDraws(end,1);
            end
            record.LastLogPost(FirstBlock)=loaded_results.logpo2(end);
            record.LastParameters(FirstBlock,:)=loaded_results.x2(end,:);
            update_last_mh_history_file(MetropolisFolder, ModelName, record);
        end
    end
    FirstBlock = find(FBlock==1,1);
end

function [d,bayestopt_,record]=set_proposal_density_to_previous_value(record,options_,bayestopt_,d,dispString)
if ~options_.use_mh_covariance_matrix
    if isfield(record,'ProposalCovariance') && isfield(record,'ProposalScaleVec')
        fprintf('%s: Recovering the previous proposal density.\n',dispString);
        d=record.ProposalCovariance;
        bayestopt_.jscale=record.ProposalScaleVec;
    else
        if ~isequal(options_.posterior_sampler_options.posterior_sampling_method,'slice')
            % pass through input d unaltered
            if options_.mode_compute~=0
                fprintf('%s: No stored previous proposal density found, continuing with the one implied by mode_compute.\n',dispString);
            elseif ~isempty(options_.mode_file)
                fprintf('%s: No stored previous proposal density found, continuing with the one implied by the mode_file.\n',dispString);
            else
                error('%s: No stored previous proposal density found, no mode-finding conducted, and no mode-file provided. I don''t know how to continue!',dispString);
            end
        end
    end
else
    % pass through input d unaltered
    fprintf('%s: ''use_mh_covariance_matrix'' specified, continuing with proposal density implied by the previous MCMC run.\n',dispString);
end

if isfield(record,'Sampler')
    if ~strcmp(record.Sampler,options_.posterior_sampler_options.posterior_sampling_method)
        warning('%s: The posterior_sampling_method %s selected differs from the %s of the original chain. This may create problems with the convergence diagnostics.',dispString,options_.posterior_sampler_options.posterior_sampling_method,record.Sampler)
        record.Sampler=options_.posterior_sampler_options.posterior_sampling_method; %update sampler used
    end
end
