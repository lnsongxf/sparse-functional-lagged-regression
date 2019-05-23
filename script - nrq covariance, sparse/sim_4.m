function sim_4(initial_seed,simulated_process,response_process,iCase) % XXX

% initial_seed=1
% simulated_process=1
% response_process=1
% iCase = 4

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



fileName_tikhonov = 'nqr_sim_tikhonov.csv';
fileName_threshold_value = 'nqr_sim_threshold_value.csv';

%% create the columns

cHeader_tikhonov = {...
    'simulated_process',...
    'simulated_process (name)',...
    'response_process',...
    'iCase',...
    'seed',...
    'exec time'...
    'CV K-fold',...
    'T',...
    'num obs',...
    'avg obs per curve',...
    'bw_mu',...
    'bw_F',...
    'bw_v',...
    'specDensity q_bartlett',...
    'specCrossDensity q_bartlett',...
    'cut_of_lag',...
    'nRandi_max',...
    'sigma^2',...
    'sigma^2 estim loc lin',...
    'sigma^2 estim loc qv (used)',...
    'cov lag 0 MSE',...
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
    'dynamic krig MSE',...
    'dynamic krig RMSE',...
    'dynamic krig MSE (true dynamics)',...
    'dynamic krig RMSE (true dynamics)',...
    'transfer: tikhonov parameter',...
    'transfer: tikhonov parameter (rule)',...
    'transfer: RMSE',...
    'Z error train',...
    'Z error train (true dynamics)',...
    'Z error test',...
    'Z error test (true dynamics)',...
    'Z error train (holdout)',...
    'Z error train (holdout true dynamics)',...
    };
cHeader_threshold_value = {...
    'simulated_process',...
    'simulated_process (name)',...
    'response_process',...
    'iCase',...
    'seed',...
    'exec time'...
    'CV K-fold',...
    'T',...
    'num obs',...
    'avg obs per curve',...
    'bw_mu',...
    'bw_F',...
    'bw_v',...
    'specDensity q_bartlett',...
    'specCrossDensity q_bartlett',...
    'cut_of_lag',...
    'nRandi_max',...
    'sigma^2',...
    'sigma^2 estim loc lin',...
    'sigma^2 estim loc qv (used)',...
    'cov lag 0 MSE',...
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
    'dynamic krig MSE',...
    'dynamic krig RMSE',...
    'dynamic krig MSE (true dynamics)',...
    'dynamic krig RMSE (true dynamics)',...
    'transfer: threshold_value',...
    'transfer: threshold_value (rule)',...
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
commaHeader_threshold_value = [cHeader_threshold_value;repmat({','},1,numel(cHeader_threshold_value))]; %insert commaas
commaHeader_threshold_value = commaHeader_threshold_value(:)';
textHeader_threshold_value = cell2mat(commaHeader_threshold_value); %cHeader in text with commas

if (~exist(fileName_tikhonov, 'file'))
    %write header to file
    fid = fopen(fileName_tikhonov,'w');
    fprintf(fid,'%s\n',textHeader_tikhonov);
    fclose(fid);
end

if (~exist(fileName_threshold_value, 'file'))
    %write header to file
    fid = fopen(fileName_threshold_value,'w');
    fprintf(fid,'%s\n',textHeader_threshold_value);
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


nGridTime = simulation_cases.nGridTime(iCase);
nRandi_max = simulation_cases.nRandi_max(iCase);


fprintf('#################################################################################\n')
fprintf('#################################################################################\n')
fprintf('\n')
fprintf('iCase = %d\n',iCase)
fprintf('seed = %d\n',initial_seed)
fprintf('nGridTime = %d\n',nGridTime)
fprintf('nRandi_max = %d\n',nRandi_max)
fprintf('simulated_process = %d\n',simulated_process)
fprintf('response_process = %d\n',response_process)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       select the tunning parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    K-FOLD CV for the common mean function

% [bw_mu, mu_est_smoother] = fKCVforMU(censor,onb);

bw_mu = 1;
mu_est_smoother = zeros(onb.nGridSmoother,1);

mu_est_ONB = onb.smoothingMatrix_gridSmoother2ONB * mu_est_smoother;
mu_est_inSpace = onb.onbMatrix*mu_est_ONB;

% subplot(1,1,1)
% plot( mu_est_inSpace )
% pause(.01)


%%    K-FOLD CV for the surface smoothing

[bw_r, bw_v, covLag0_smoother, covLag0_ridge_smoother, sigma2_est_loc_qv, sigma2_est_loc_lin] = fKCVforR_locQv(censor,onb,mu_est_inSpace);

covLag0_ONB = onb.smoothingMatrix_gridSmoother2ONB * covLag0_smoother * onb.smoothingMatrix_gridSmoother2ONB';
covLag0_inSpace = onb.onbMatrix * covLag0_ONB * onb.onbMatrix';

% subplot(1,1,1)
% surf( covLag0_inSpace )
% pause(.01)

%%  K-FOLD CV for the cross-spectral density
[bw_c, crossCov0_smoother] = fKCVforCrossR(censor,onb,response,mu_est_inSpace);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           estimate the spectral density and quantify the estimation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% estimate the spectral density
bw_F = bw_r; % the same smoothing parameter as for cov lag 0
specDensity = fEstimateSpecDensity_v2(censor,onb,mu_est_inSpace,bw_F,nGridFreq);

%% extract auto-covariances (integrate out the spectral density)
numOfLags_all = censor.nGridTime - 1;
covLags = fEstimateCovLags(onb,specDensity,numOfLags_all);

%% quantify the error of spec density estimation
specDensity_error = fSpecDensityError( onb, specDensity, specDensity_true );
% specDensity_error.mse_fro
% specDensity_error.rmse_fro
% specDensity_error.r2mse_fro
% specDensity_error.mse_intmax
% specDensity_error.rmse_intmax
% specDensity_error.r2mse_intmax
% specDensity_error.mse_maxmax
% specDensity_error.rmse_maxmax

%% quantify the error of cov lag 0 estimation
covLag0_error = fCovLag0Error( onb, covsLags_true, covLag0_ONB );
% covLag0_error.mse
% covLag0_error.rmse

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
% %%
% pause(0.01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     kriging for the process X

% if I want to get confidence bands
cofidence_bands = 0;

% use the sigma from the loc-qv estimation
s2 = sigma2_est_loc_qv;
if s2 < 1e-2, s2 = 1e-2; end
sigma_est = sqrt( s2 );


% kriging on the TRAINING sample
disp('kriging 1./4')
[krigingX_dyna,krigingX_dyna_padded] = fKrigingX( censor, onb, covLags, mu_est_ONB, sigma_est, 1, cofidence_bands, dataFull_padding );
disp('kriging 2./4')
[krigingX_dyna_true,krigingX_dyna_padded_true] = fKrigingX( censor, onb, covsLags_true, mu_true, sigma, 1, cofidence_bands, dataFull_padding );

% kriging on the TEST sample
disp('kriging 3./4')
[krigingX_dyna_test,krigingX_dyna_test_padded] = fKrigingX( censor_test, onb, covLags, mu_est_ONB, sigma_est, 1, cofidence_bands, dataFull_padding );
disp('kriging 4./4')
[krigingX_dyna_test_true,krigingX_dyna_test_padded_true] = fKrigingX( censor_test, onb, covsLags_true, mu_true, sigma, 1, cofidence_bands, dataFull_padding );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   quantify the kriging error

krigingX_dyna_error = fKrigingX_error( onb, dataFull, krigingX_dyna, covsLags_true.lag0, cofidence_bands );
krigingX_dyna_error_true = fKrigingX_error( onb, dataFull, krigingX_dyna_true, covsLags_true.lag0, cofidence_bands );
% krigingX_stat_error = fKrigingX_error( onb, dataFull, krigingX_stat, covsLags_true.lag0, cofidence_bands );
% krigingX_stat_error_true = fKrigingX_error( onb, dataFull, krigingX_stat_true, covsLags_true.lag0, cofidence_bands );

% xxx.mse
% xxx.rmse
% xxx.coverage_point
% xxx.coverage_simul

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           estimate the CROSS-spectral density
disp('cross spectral density')
response_mean = 0;
specCrossDensity = fEstimateSpecCrossDensity_holdout_v5(censor,onb,response,response_mean,mu_est_inSpace,bw_c,nGridFreq,holdout); % bandwidth: bw_c

specCrossDensity_error = fSpecCrossDensityError( onb, specCrossDensity, specCrossDensity_true );
% specCrossDensity_error.mse_qv
% specCrossDensity_error.rmse_qv
% specCrossDensity_error.mse_intmax
% specCrossDensity_error.rmse_intmax


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the spectral transfer function

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
% 
% 
% %%
% xxx
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tikhonov regularization

tikhonovParameter = 2* censor.nGridTime ^(-1/3) * (sum(censor.nGrid)/nGridTime)^(-1/3);
tikhonovParameter_string = '2* censor.nGridTime ^(-1/3) * (sum(censor.nGrid)/nGridTime)^(-1/3)';
disp( ['spectral transfer, tikhonov parameter = ',num2str(tikhonovParameter)] )
transfer = fEstimateTransfer_Tikhonov(onb,specDensity,specCrossDensity, tikhonovParameter);

%% quantify the error of transfer functionals
rmse_transfer = fTransferError(transfer,transfer_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the response process

% estimation on the TRAINING sample
response_est = fKrigingZ(krigingX_dyna_padded,transfer,dataFull_padding);
response_est_true = fKrigingZ(krigingX_dyna_padded_true,transfer_true,dataFull_padding);

% estimation on the TEST sample
response_est_test = fKrigingZ(krigingX_dyna_test_padded,transfer,dataFull_padding);
response_est_test_true = fKrigingZ(krigingX_dyna_test_padded_true,transfer_true,dataFull_padding);

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exectime = toc;

results = {...
    num2str( simulated_process ),...
    'based on quadratic rational',...
    num2str( response_process ),...
    num2str( iCase ),...
    num2str( initial_seed ),...
    num2str( exectime ),...
    num2str( censor.cv_K ),...
    num2str( censor.nGridTime ),...
    num2str( sum(censor.nGrid) ),...
    num2str( sum(censor.nGrid)/nGridTime ),...
    num2str( bw_mu ),...
    num2str( bw_F ),...
    num2str( bw_v ),...
    num2str( specDensity.q_bartlett ),...
    num2str( specCrossDensity.q_bartlett ),...
    num2str( covLags.cutOfLag ),...
    num2str( nRandi_max ),...
    num2str( sigma^2 ),...
    num2str( sigma2_est_loc_lin ),...
    num2str( sigma2_est_loc_qv ),...
    num2str( covLag0_error.mse ),...
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
    num2str( krigingX_dyna_error.mse ),...
    num2str( krigingX_dyna_error.rmse ),...
    num2str( krigingX_dyna_error_true.mse ),...
    num2str( krigingX_dyna_error_true.rmse ),...
    num2str( tikhonovParameter ),...
    num2str( tikhonovParameter_string ),...    
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

threshold_value = censor.nGridTime.^(-1/3);
threshold_value_string = 'censor.nGridTime.^(-1/3)';
disp( ['spectral transfer, threshold_value= ',num2str(threshold_value)] )
transfer = fEstimateTransfer_threshold_value(onb,specDensity,specCrossDensity, threshold_value);

%% quantify the error of transfer functionals
rmse_transfer = fTransferError(transfer,transfer_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the response process

% estimation on the TRAINING sample
response_est = fKrigingZ(krigingX_dyna_padded,transfer,dataFull_padding);
response_est_true = fKrigingZ(krigingX_dyna_padded_true,transfer_true,dataFull_padding);

% estimation on the TEST sample
response_est_test = fKrigingZ(krigingX_dyna_test_padded,transfer,dataFull_padding);
response_est_test_true = fKrigingZ(krigingX_dyna_test_padded_true,transfer_true,dataFull_padding);

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
    'based on quadratic rational',...
    num2str( response_process ),...
    num2str( iCase ),...
    num2str( initial_seed ),...
    num2str( exectime ),...
    num2str( censor.cv_K ),...
    num2str( censor.nGridTime ),...
    num2str( sum(censor.nGrid) ),...
    num2str( sum(censor.nGrid)/nGridTime ),...
    num2str( bw_mu ),...
    num2str( bw_F ),...
    num2str( bw_v ),...
    num2str( specDensity.q_bartlett ),...
    num2str( specCrossDensity.q_bartlett ),...
    num2str( covLags.cutOfLag ),...
    num2str( nRandi_max ),...
    num2str( sigma^2 ),...
    num2str( sigma2_est_loc_lin ),...
    num2str( sigma2_est_loc_qv ),...
    num2str( covLag0_error.mse ),...
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
    num2str( krigingX_dyna_error.mse ),...
    num2str( krigingX_dyna_error.rmse ),...
    num2str( krigingX_dyna_error_true.mse ),...
    num2str( krigingX_dyna_error_true.rmse ),...
    num2str( threshold_value ),...
    num2str( threshold_value_string ),...
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
fid = fopen(fileName_threshold_value,'a');
fprintf(fid,'%s\n',textResults);
fclose(fid);


end % end of function
