options_.debug=true;

options_.lik_init=1;
estimation(order=1,datafile='../fsdat_simul',nobs=192,mode_file=@{mode_file_name},mode_compute=0,kalman_algo=1,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous);
fval_algo_1=oo_.likelihood_at_initial_parameters;
if isfield(oo_,'SmoothedMeasurementErrors')
    SmoothedMeasurementErrors(:,:,1)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
end
SmoothedShocks(:,:,1)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,1)=cell2mat(struct2cell(oo_.SmoothedVariables));

options_.lik_init=1;
estimation(order=1,datafile='../fsdat_simul',nobs=192,mode_file=@{mode_file_name},mode_compute=0,kalman_algo=2,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous);
fval_algo_2=oo_.likelihood_at_initial_parameters;
if isfield(oo_,'SmoothedMeasurementErrors')
    SmoothedMeasurementErrors(:,:,2)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
end
SmoothedShocks(:,:,2)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,2)=cell2mat(struct2cell(oo_.SmoothedVariables));

options_.lik_init=1;
estimation(order=1,datafile='../fsdat_simul',nobs=192,mode_file=@{mode_file_name},mode_compute=0,kalman_algo=3,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous);
fval_algo_3=oo_.likelihood_at_initial_parameters;
if isfield(oo_,'SmoothedMeasurementErrors')
    SmoothedMeasurementErrors(:,:,3)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
end
SmoothedShocks(:,:,3)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,3)=cell2mat(struct2cell(oo_.SmoothedVariables));

options_.lik_init=1;
estimation(order=1,datafile='../fsdat_simul',nobs=192,mode_file=@{mode_file_name},mode_compute=0,kalman_algo=4,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous);
fval_algo_4=oo_.likelihood_at_initial_parameters;
if isfield(oo_,'SmoothedMeasurementErrors')
    SmoothedMeasurementErrors(:,:,4)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
end
SmoothedShocks(:,:,4)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,4)=cell2mat(struct2cell(oo_.SmoothedVariables));

if isfield(oo_,'SmoothedMeasurementErrors')
    if max(max(abs(SmoothedMeasurementErrors-repmat(SmoothedMeasurementErrors(:,:,1),[1,1,4]))))>1e-6
        error('SmoothedMeasurementErrors do not match')
    end
end

if max(max(abs(SmoothedShocks-repmat(SmoothedShocks(:,:,1),[1,1,4]))))>1e-6
    error('SmoothedShocks do not match')
end

if max(max(abs(SmoothedVariables-repmat(SmoothedVariables(:,:,1),[1,1,4]))))>1e-6
    error('SmoothedVariables do not match')
end

if max(abs([fval_algo_2,fval_algo_3,fval_algo_4]-fval_algo_1))>1e-5
    error('Likelihoods do not match')
end