function sim_3_full(initial_seed,simulated_process,response_process,nGridTime) % XXX

% initial_seed=1
% simulated_process=4
% response_process=1
% nGridTime = 300

nRandi_max = 10; % this is only to use the existing fSimulate_XX function


% simulated_process: 0=iid, 1=AR(1,0.7), 4=MA(4)
% response_process:   1=(b_0,b_1),  2=(b_0,b_3),  3=slowly decay   all with sigma^2 = 0.1
% iCase: 1..15

holdout = 0.2; % set up the holdout parameter, the holdout partition is not used for fitting the cross-spectral density, it's used for error evaluation

prepare_simulation_cases;

%% random seeds
rng(initial_seed) % XXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the CSV file for saving the results



% fileName = ['z_sim_01_c', num2str(iCase,'%02d'),'_p',num2str(simulated_process),'_r',num2str(response_process,'%02d'), '.csv'];
fileName_tikhonov = 'nrq_sim_complete_tikhonov.csv';
fileName_threshold_multiplier = 'nrq_sim_complete_threshold.csv';

%% create the columns

cHeader_tikhonov = {...
    'simulated_process',...
    'response_process',...    
    'seed',...
    'exec time'...    
    'T',...    
    'q_bartlett',...    
    'cov lag 0 RMSE',...
    'spec density MSE - FRO',...
    'spec density RMSE - FRO',...
    'spec density R2MSE - FRO',...
    'spec density MSE - intmax',...
    'spec density RMSE - intmax',...
    'spec density R2MSE - intmax',...
    'spec density MSE - maxmax',...
    'spec density RMSE - maxmax',...
    'spec cross density MSE - L2',...
    'spec cross density RMSE - L2',...
    'spec cross density MSE - intmax',...
    'spec cross density RMSE - intmax',...    
    'transfer: tikhonov parameter',...
    'transfer: RMSE',...
    'Z error train',...
    'Z error train (true dynamics)',...
    'Z error test',...
    'Z error test (true dynamics)',...
    'Z error train (holdout)',...
    'Z error train (holdout true dynamics)',...
    };
cHeader_threshold_multiplier = {...
    'simulated_process',...
    'response_process',...    
    'seed',...
    'exec time'...    
    'T',...    
    'q_bartlett',...    
    'cov lag 0 RMSE',...
    'spec density MSE - FRO',...
    'spec density RMSE - FRO',...
    'spec density R2MSE - FRO',...
    'spec density MSE - intmax',...
    'spec density RMSE - intmax',...
    'spec density R2MSE - intmax',...
    'spec density MSE - maxmax',...
    'spec density RMSE - maxmax',...
    'spec cross density MSE - L2',...
    'spec cross density RMSE - L2',...
    'spec cross density MSE - intmax',...
    'spec cross density RMSE - intmax',...    
    'transfer: truncation parameter',...
    'transfer: RMSE',...
    'Z error train',...
    'Z error train (true dynamics)',...
    'Z error test',...
    'Z error test (true dynamics)',...
    'Z error train (holdout)',...
    'Z error train (holdout true dynamics)',...
    };


%%
commaHeader_tikhonov = [cHeader_tikhonov;repmat({','},1,numel(cHeader_tikhonov))]; %insert commaas
commaHeader_tikhonov = commaHeader_tikhonov(:)';
textHeader_tikhonov = cell2mat(commaHeader_tikhonov); %cHeader in text with commas

commaHeader_threshold_multiplier = [cHeader_threshold_multiplier;repmat({','},1,numel(cHeader_threshold_multiplier))]; %insert commaas
commaHeader_threshold_multiplier = commaHeader_threshold_multiplier(:)';
textHeader_threshold_multiplier = cell2mat(commaHeader_threshold_multiplier); %cHeader in text with commas

if (~exist(fileName_tikhonov, 'file'))
    %write header to file
    fid = fopen(fileName_tikhonov,'w');
    fprintf(fid,'%s\n',textHeader_tikhonov);
    fclose(fid);
end

if (~exist(fileName_threshold_multiplier, 'file'))
    %write header to file
    fid = fopen(fileName_threshold_multiplier,'w');
    fprintf(fid,'%s\n',textHeader_threshold_multiplier);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the onb.basis

% nGridSpace = 241, nKnots = 24
onb = fPrepareONB(241, 24);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                          SIMULATION STARTS HERE                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% random seed that has been saved
tic



fprintf('#################################################################################\n')
fprintf('#################################################################################\n')
fprintf('\n')
fprintf('seed = %d\n',initial_seed)
fprintf('simulated_process = %d\n',simulated_process)
fprintf('response_process = %d\n',response_process)
fprintf('nGridTime = %d\n',nGridTime)
fprintf('\n')
fprintf('#################################################################################\n')
fprintf('#################################################################################\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate the process - twice
nGridFreq = 1000;
dataFull_padding = 10;

% if simulated_process <= 2
%     % AR(1)
%     pars_simulated_process = fSimulate_AR_random_prepare_coef(onb,simulated_process);
%     [censor, dataFull,dataFull_padded, covsLags_true, mu_true, specDensity_true, maxEigenvalue,sigma] = fSimulate_AR_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
%     [censor_test, dataFull_test,dataFull_test_padded,~,~,~,~,~] = fSimulate_AR_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
% else
%     % MA(q)
%     pars_simulated_process = fSimulate_MA_random_prepare_coef(onb);
%     [censor, dataFull,dataFull_padded, covsLags_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
%     [censor_test, dataFull_test,dataFull_test_padded,~,~,~,~,~] = fSimulate_MA_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
% end


if simulated_process <= 2
    % AR(1)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, maxEigenvalue,sigma] = fSimulate_AR(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_AR(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
else
    % MA(q)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_MA(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
end



%% set up the response process
pars_response_process = fSimulate_response_prepare_coef(onb.nGridSpace,10); % prepare 10 random functions
[response, transfer_true, specCrossDensity_true] = fSimulate_response_ONB(onb,response_process,pars_response_process,dataFullONB_padded,dataFull_padding,specDensity_true);
[response_test,~,~] = fSimulate_response_ONB(onb,response_process,pars_response_process,dataFullONB_test_padded,dataFull_padding,specDensity_true);

%% visualization

% for t = 1:12
%     subplot(3,4,t)
%     plot(onb.gridSpace, dataFull(t,:))
%     hold on
%     plot( onb.gridSpace(censor.grid{t}) , censor.data{t},'x' )
%     hold off
%     integral = mean(dataFull(t,:));
%     ylim([-4 4])
%     xlim([0 1])
%     title( ['t=',num2str(t),', int=',num2str(round(integral,2)),', Y=',num2str(round(response(t),2))] )
% end
% suptitle('first 12 snapshots')
% pause(0.1)


%
% %%
% subplot(2,1,1)
% m = onb.onbMatrix *transfer_true.specTransfer';
% plot( real( m(1,:) ) )
% title('real')
% subplot(2,1,2)
% plot( imag( m(1,:) ) )
% title('imag')
%
% xxx
% %%
% subplot(2,1,1)
% surf( real( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )
% subplot(2,1,2)
% surf( imag( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% working with variables:
%  - dataFull .... all curves, dataFull(1,:) is the first curve in spatial variable
%  - dataFullONB .... all curves, dataFull(1,:) is the first curve in spatial variable
%  - response .... response process




% mind the fact that the crosscov is done only without holdout
% holdout = percentage of how much data at the end of the response process is not used for the training
% typically, holdout = 0.8
holdout_number = floor(censor.nGridTime * holdout); % the number of points in the holdout partition
holdout_start = censor.nGridTime - holdout_number + 1; % the index of the first datum in the holdout partition

% prepare the estimates of the autocovariance operators
num_of_lags = nGridTime;
cov_estim = nan(num_of_lags,onb.nBasis,onb.nBasis);
crosscov_estim = nan(num_of_lags,onb.nBasis);
crosscov_estim_minus_lags = nan(num_of_lags,onb.nBasis);
for lag_real = 0:(num_of_lags-1)
    lag_index = lag_real + 1;
    
    cov_estim(lag_index,:,:) = 1/(nGridTime - lag_real) * dataFullONB( (lag_real+1):nGridTime ,:)' * dataFullONB( 1:(nGridTime-lag_real) ,:);
    crosscov_estim(lag_index,:) = 1/((holdout_start-1) - lag_real) * response( (lag_real+1):(holdout_start-1))' * dataFullONB( 1:((holdout_start-1)-lag_real) ,:);
    crosscov_estim_minus_lags(lag_index,:) = 1/((holdout_start-1) - lag_real) * response( 1:((holdout_start-1)-lag_real) )' * dataFullONB( (lag_real+1):(holdout_start-1) ,:);
end



%% quantify the error of cov lag 0 estimation
covLag0_error = fCovLag0Error( onb, covsLags_true, squeeze(cov_estim(1,:,:)) );
%%

% surf( onb.onbMatrix * squeeze(cov_estim(60,:,:)) * onb.onbMatrix'  )
% surf( onb.onbMatrix * squeeze(covsLags_true.ONB(5,:,:)) * onb.onbMatrix' )
% surf( onb.onbMatrix * crosscov_estim(1:100,:)' )

%% estimate the spectral density
nGridFreq = 1000;
specDensity = [];
specDensity.nGridFreq = nGridFreq;
specDensity.gridFreq = linspace(-pi,pi,specDensity.nGridFreq);
q_bartlett = floor(sqrt(nGridTime));
specDensity.q_bartlett = q_bartlett;
specDensity.ONB = nan(nGridFreq, onb.nBasis, onb.nBasis);

specCrossDensity = [];
specCrossDensity.nGridFreq = nGridFreq;
specCrossDensity.gridFreq = linspace(-pi,pi, specCrossDensity.nGridFreq);
specCrossDensity.ONB = nan(specCrossDensity.nGridFreq, onb.nBasis);
specCrossDensity.q_bartlett = q_bartlett;

for k = 1:nGridFreq    
    omega = specDensity.gridFreq(k);
    % zero lag
    specDensity.ONB(k,:,:) = 1/(2*pi)*squeeze(cov_estim(1,:,:));
    specCrossDensity.ONB(k,:) = 1/(2*pi)*crosscov_estim(1,:);
    % non-zero lags
    for lag_real = 1:(q_bartlett-1)
        lag_index = lag_real + 1;
        specDensity.ONB(k,:,:) = squeeze(specDensity.ONB(k,:,:)) +... % positive lag
            1/(2*pi)*(1 - lag_real/q_bartlett)*squeeze(cov_estim(lag_index,:,:)) * exp(-1i*omega*lag_real);
        specDensity.ONB(k,:,:) = squeeze(specDensity.ONB(k,:,:)) +... % negative lag => transpose the autocovariance
            1/(2*pi)*(1 - lag_real/q_bartlett)*squeeze(cov_estim(lag_index,:,:))' * exp(+1i*omega*lag_real);
        specCrossDensity.ONB(k,:) = specCrossDensity.ONB(k,:) +...
            1/(2*pi)*(1 - lag_real/q_bartlett)*crosscov_estim(lag_index,:) * exp(-1i*omega*lag_real);
        specCrossDensity.ONB(k,:) = specCrossDensity.ONB(k,:) +...
            1/(2*pi)*(1 - lag_real/q_bartlett)*crosscov_estim_minus_lags(lag_index,:) * exp(+1i*omega*lag_real);
    end
    
end


%% quantify the error of spec density estimation
specDensity_error = fSpecDensityError( onb, specDensity, specDensity_true );

 
% 
% %% visualize spectral density estimation
% traces = nan(1,specDensity.nGridFreq);
% ranks = nan(1,specDensity.nGridFreq);
% traces_true = nan(1,specDensity.nGridFreq);
% ranks_true = nan(1,specDensity.nGridFreq);
% for k = 1:specDensity.nGridFreq
%     m1 = squeeze( specDensity.ONB(k,:,:) );
%     m2 = squeeze( specDensity_true.ONB(k,:,:) );
%     traces(k) = trace(m1);
%     traces_true(k) = trace(m2);
%     
%     e1 = eig(m1);
%     ranks(k) = sum( e1 > sum(e1)*0.01 );
%     e2 = eig(m2);
%     ranks_true(k) = sum( e2 > sum(e2)*0.01 );
% end
% 
% subplot(2,1,1)
% plot(traces)
% hold on
% plot(real(traces_true))
% hold off
% legend('estim','true')
% 
% subplot(2,1,2)
% plot(ranks)
% hold on
% plot(ranks_true)
% hold off
% legend('estim','true')
% suptitle('spectral density estimation')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           estimate the CROSS-spectral density
disp('cross spectral density')
response_mean = 0;
% specCrossDensity = fEstimateSpecCrossDensity_holdout(censor,onb,response,response_mean,mu_est_inSpace,bw_c,nGridFreq,holdout); % bandwidth: bw_c

specCrossDensity_error = fSpecCrossDensityError( onb, specCrossDensity, specCrossDensity_true );
% specCrossDensity_error.mse_qv
% specCrossDensity_error.rmse_qv
% specCrossDensity_error.mse_intmax
% specCrossDensity_error.rmse_intmax

% 
% % vizualization of specCrossDensity
% subplot(2,2,1)
% surf( transfer_true.gridFreq, onb.gridSpace, real(onb.onbMatrix * specCrossDensity_true.ONB'), 'EdgeColor','none' )
% title('REAL part, true')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('real part')
% subplot(2,2,2)
% surf( transfer_true.gridFreq, onb.gridSpace, imag(onb.onbMatrix * specCrossDensity_true.ONB'), 'EdgeColor','none' )
% title('IMAG part, true')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('imag part')
% 
% subplot(2,2,3)
% surf( transfer_true.gridFreq, onb.gridSpace, real(onb.onbMatrix * specCrossDensity.ONB'), 'EdgeColor','none' )
% title('REAL part, estim')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('real part')
% subplot(2,2,4)
% surf( transfer_true.gridFreq, onb.gridSpace, imag(onb.onbMatrix * specCrossDensity.ONB'), 'EdgeColor','none' )
% title('IMAG part, estim')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('imag part')
% suptitle('spectral cross-density')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % estimate the spectral transfer function
% 
% tikhonovParameter=0.15;
% disp( ['spectral transfer, tikhonov parameter = ',num2str(tikhonovParameter)] )
% transfer = fEstimateTransfer_Tikhonov(onb,specDensity,specCrossDensity, tikhonovParameter);
% % transfer = fEstimateTransfer_Tichonov(onb,specDensity,specCrossDensity, 0.1);
% % transfer = fEstimateTransfer(onb,specDensity_true,specCrossDensity_true, truncationEigen);
% 
% % vizualization of spectral transfer functions
% subplot(2,2,1)
% surf( transfer_true.gridFreq, onb.gridSpace, real( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )
% title('REAL part, true')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('real part')
% subplot(2,2,2)
% surf( transfer_true.gridFreq, onb.gridSpace, imag( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )
% title('IMAG part, true')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('imag part')
% 
% subplot(2,2,3)
% surf( transfer.gridFreq, onb.gridSpace, real( onb.onbMatrix *transfer.specTransfer' ), 'EdgeColor','none' )
% title('REAL part, estim')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('real part')
% subplot(2,2,4)
% surf( transfer.gridFreq, onb.gridSpace, imag( onb.onbMatrix *transfer.specTransfer' ), 'EdgeColor','none' )
% title('IMAG part, estim')
% xlabel('frequency')
% ylabel('space [0,1]')
% zlabel('imag part')
% 
% suptitle('spectral transfer function')
% 
% %% visualize the estimated regression functionals
% for lag_real = -5:5
%     lag_index = lag_real + transfer.numOfLags+1;
%     subplot(3,4,lag_real+6)
%     plot( onb.onbMatrix* transfer.ONB(lag_index,:)' )
%     hold on
%     plot( onb.onbMatrix* transfer_true.ONB(lag_index,:)')
%     hold off
%     title( ['lag = ',num2str(lag_real)])
%     legend('estim','true')
%     ylim([-2,2])
% end
% % 
% % 
% 



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Tikhonov regularization

tikhonovParameter = 1 / specDensity.q_bartlett;
disp( ['spectral transfer, tikhonov parameter = ',num2str(tikhonovParameter)] )
transfer = fEstimateTransfer_Tikhonov(onb,specDensity,specCrossDensity, tikhonovParameter);

%% quantify the error of transfer functionals
rmse_transfer = fTransferError(transfer,transfer_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the response process


% here I'm especially padding with zeroes because that's what Hoermann does
dataFullONB_padded_w_zeros=[];
dataFullONB_padded_w_zeros.est = [zeros(dataFull_padding,onb.nBasis); dataFullONB; zeros(dataFull_padding,onb.nBasis)];
dataFullONB_test_padded_w_zeros=[];
dataFullONB_test_padded_w_zeros.est = [zeros(dataFull_padding,onb.nBasis); dataFullONB_test; zeros(dataFull_padding,onb.nBasis)];

% estimation on the TRAINING sample
response_est = fKrigingZ(dataFullONB_padded_w_zeros,transfer,dataFull_padding);
response_est_true = fKrigingZ(dataFullONB_padded_w_zeros,transfer_true,dataFull_padding);

% estimation on the TEST sample
response_est_test = fKrigingZ(dataFullONB_test_padded_w_zeros,transfer,dataFull_padding);
response_est_test_true = fKrigingZ(dataFullONB_test_padded_w_zeros,transfer_true,dataFull_padding);

% calculate var(Z) for normalization
var_z = 0;
for kk = 1:(2*transfer_true.numOfLags-1)
    for ll = 1:(2*transfer_true.numOfLags-1)
        lag = kk-ll;
        if lag >= 0
            C = squeeze(covsLags_true.ONB(lag+1,:,:));
        else
            C = squeeze(covsLags_true.ONB(-lag+1,:,:))';
        end
        var_z = var_z + transfer_true.ONB(kk,:)*C*transfer_true.ONB(ll,:)';
    end
end

% quantify errors
Z_est_error = fKrigingZ_error( response_est, response', var_z );
Z_est_true_error = fKrigingZ_error( response_est_true, response', var_z );
Z_est_test_error = fKrigingZ_error( response_est_test, response_test', var_z );
Z_est_test_true_error = fKrigingZ_error( response_est_test_true, response_test', var_z );
Z_est_holdout_error = fKrigingZ_holdout_error( response_est, response', var_z, holdout);
Z_est_holdout_true_error = fKrigingZ_holdout_error( response_est_true, response', var_z, holdout );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exectime = toc;

results = {...
    num2str( simulated_process ),...
    num2str( response_process ),...
    num2str( initial_seed ),...
    num2str( exectime ),...
    num2str( censor.nGridTime ),...
    num2str( specDensity.q_bartlett ),...
    num2str( covLag0_error.rmse ),...
    num2str( specDensity_error.mse_fro ),...
    num2str( specDensity_error.rmse_fro ),...
    num2str( specDensity_error.r2mse_fro ),...
    num2str( specDensity_error.mse_intmax ),...
    num2str( specDensity_error.rmse_intmax ),...
    num2str( specDensity_error.r2mse_intmax ),...
    num2str( specDensity_error.mse_maxmax ),...
    num2str( specDensity_error.rmse_maxmax ),...
    num2str( specCrossDensity_error.mse_qv ),...
    num2str( specCrossDensity_error.rmse_qv ),...
    num2str( specCrossDensity_error.mse_intmax ),...
    num2str( specCrossDensity_error.rmse_intmax ),...
    num2str( tikhonovParameter ),...
    num2str( rmse_transfer ),...
    num2str( Z_est_error ),...
    num2str( Z_est_true_error ),...
    num2str( Z_est_test_error ),...
    num2str( Z_est_test_true_error ),...
    num2str( Z_est_holdout_error ),...
    num2str( Z_est_holdout_true_error )...
    };

% write the results into CSV
commaResults = [results;repmat({','},1,numel(results))]; %insert commaas
commaResults = commaResults(:)';
textResults = cell2mat(commaResults); %cHeader in text with commas

% write into seperate file for each triplet: iCase, simulated_process, response_process
%         fid = fopen(fileName,'a');
%         fprintf(fid,'%s\n',textResults);
%         fclose(fid);

% write into one file
fid = fopen(fileName_tikhonov,'a');
fprintf(fid,'%s\n',textResults);
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Truncation regularization - varying over frequencies

threshold_multiplier = 1;
disp( ['spectral transfer, threshold_multiplier= ',num2str(threshold_multiplier)] )
transfer = fEstimateTransfer_threshold(onb,specDensity,specCrossDensity, threshold_multiplier);

%% quantify the error of transfer functionals
rmse_transfer = fTransferError(transfer,transfer_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the response process

% estimation on the TRAINING sample
response_est = fKrigingZ(dataFullONB_padded_w_zeros,transfer,dataFull_padding);
response_est_true = fKrigingZ(dataFullONB_padded_w_zeros,transfer_true,dataFull_padding);

% estimation on the TEST sample
response_est_test = fKrigingZ(dataFullONB_test_padded_w_zeros,transfer,dataFull_padding);
response_est_test_true = fKrigingZ(dataFullONB_test_padded_w_zeros,transfer_true,dataFull_padding);

% calculate var(Z) for normalization
var_z = 0;
for kk = 1:(2*transfer_true.numOfLags-1)
    for ll = 1:(2*transfer_true.numOfLags-1)
        lag = kk-ll;
        if lag >= 0, C = squeeze(covsLags_true.ONB(lag+1,:,:));
        else, C = squeeze(covsLags_true.ONB(-lag+1,:,:))'; end
        var_z = var_z + transfer_true.ONB(kk,:)*C*transfer_true.ONB(ll,:)';
    end
end

% quantify errors
Z_est_error = fKrigingZ_error( response_est, response', var_z );
Z_est_true_error = fKrigingZ_error( response_est_true, response', var_z );
Z_est_test_error = fKrigingZ_error( response_est_test, response_test', var_z );
Z_est_test_true_error = fKrigingZ_error( response_est_test_true, response_test', var_z );
Z_est_holdout_error = fKrigingZ_holdout_error( response_est, response', var_z, holdout);
Z_est_holdout_true_error = fKrigingZ_holdout_error( response_est_true, response', var_z, holdout );




%% save results

exectime = toc;


results = {...
    num2str( simulated_process ),...
    num2str( response_process ),...
    num2str( initial_seed ),...
    num2str( exectime ),...
    num2str( censor.nGridTime ),...
    num2str( specDensity.q_bartlett ),...
    num2str( covLag0_error.rmse ),...
    num2str( specDensity_error.mse_fro ),...
    num2str( specDensity_error.rmse_fro ),...
    num2str( specDensity_error.r2mse_fro ),...
    num2str( specDensity_error.mse_intmax ),...
    num2str( specDensity_error.rmse_intmax ),...
    num2str( specDensity_error.r2mse_intmax ),...
    num2str( specDensity_error.mse_maxmax ),...
    num2str( specDensity_error.rmse_maxmax ),...
    num2str( specCrossDensity_error.mse_qv ),...
    num2str( specCrossDensity_error.rmse_qv ),...
    num2str( specCrossDensity_error.mse_intmax ),...
    num2str( specCrossDensity_error.rmse_intmax ),...
    num2str( threshold_multiplier ),...
    num2str( rmse_transfer ),...
    num2str( Z_est_error ),...
    num2str( Z_est_true_error ),...
    num2str( Z_est_test_error ),...
    num2str( Z_est_test_true_error ),...
    num2str( Z_est_holdout_error ),...
    num2str( Z_est_holdout_true_error )...
    };

% write the results into CSV
commaResults = [results;repmat({','},1,numel(results))]; %insert commaas
commaResults = commaResults(:)';
textResults = cell2mat(commaResults); %cHeader in text with commas

% write into seperate file for each triplet: iCase, simulated_process, response_process
%         fid = fopen(fileName,'a');
%         fprintf(fid,'%s\n',textResults);
%         fclose(fid);

% write into one file
fid = fopen(fileName_threshold_multiplier,'a');
fprintf(fid,'%s\n',textResults);
fclose(fid);

