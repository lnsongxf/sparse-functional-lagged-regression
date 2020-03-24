%{
This is a DEMO program demonstrating the functional lagged regression
methodology from sparse and noisy observations. The program:
 1) simulates the data from a specified process and regression scheme,
 2) estimates the dynamics of the process,
 3) displays predictions of the latent functional regressor data,
 4) estimates the filter coefficients by truncation and Tikhonov
 regularisation,
 5) displays predictions of the response process.

Visualisations are plotted using MATLAB's plot functions and the program
displays the figures while calculating.

Principal references:
   Rubin, T. & Panaretos, V. M. (2019). Functional Lagged Regression with Sparse Noisy Observations. arXiv:1905.07218.

Link for GitHub repository via DOI:
   http://doi.org/10.5281/zenodo.3190741

Code developed on: MATLAB R2018a

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    1)      set up the parameters of the simulated processes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% random seed for the RNG
initial_seed = 1234;

%%% set up the dynamics of the regressor FTS {X_t}
simulated_process = 1; % 1 = FAR(1)
% simulated_process = 4; % 2 = FMA(4)

%%% set up the filter coefficients shape
% response_process_shape = 1; % 1 = constant one function
% response_process_shape = 2; % 2 = cos( 4*pi*x ) ... option (A) from the paper
response_process_shape = 3; % 3 = sin( 2*pi*x ) ..... option (B) from the paper
% response_process_shape = 4 % 4 = ( -x+0.5)*exp( -20*( -x+0.5)^2 ) 

%%% set up the regression structure
% response_process = 1; % 1 = (b_0,b_1) ............ (reg1) scheme from the paper
% response_process = 2; % 2 = (b_0,b_3) .......... (reg2) scheme from the paper
response_process = 3; % 3 = slowly decaying .... (reg3) scheme from the paper

%%% sample size parameter - read the WARNING below
% iCase =  1; % T =  300, N^max = 10
% iCase =  2; % T =  600, N^max = 10
% iCase =  3; % T =  900, N^max = 10
% iCase =  4; % T = 1200, N^max = 10
% iCase =  5; % T =  300, N^max = 20
iCase =  6; % T =  600, N^max = 20
% iCase =  7; % T =  900, N^max = 20
% iCase =  8; % T = 1200, N^max = 20
% iCase =  9; % T =  300, N^max = 40
% iCase = 10; % T =  600, N^max = 40
% iCase = 11; % T =  900, N^max = 40
% iCase = 12; % T = 1200, N^max = 40
% iCase = 13; % T =  300, N^max = 60
% iCase = 14; % T =  600, N^max = 60
% iCase = 15; % T =  900, N^max = 60
% iCase = 16; % T = 1200, N^max = 60

% WARNING: sample sizes with T >= 900 or N^max >= 40 take long to run and
% the runtime needs lots of RAM. There is a lot of potential to optimize
% the code for both speed and memory required.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    2)      prepare simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the random seed was set by user above
rng(initial_seed)

% here's the library
addpath('master')
addpath('master/fdaM')

% set up the holdout parameter, the holdout partition is not used for fitting the cross-spectral density, it's used for error evaluation
holdout = 0.2;

% prepare the ONB basis
onb = fPrepareONB(241, 24);

% access the sample size parameters, given "iCase" set up by the user above
[simulation_cases,~,~,~,~] = fPrepareSimulationCases();
nRandi_min = 0;
nGridTime = simulation_cases.nGridTime(iCase);
nRandi_max = simulation_cases.nRandi_max(iCase);

% the resolution of the frequency grid on [-pi,pi]
nGridFreq = 1000;

% the padding of the simulated process before and after 1...T for plugging
% into the filter coefficients
dataFull_padding = 30;

% simulate the process - once for the actual simulation, the other for
% prediction error testing on an independet sample
if simulated_process <= 2
    % AR(1)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, maxEigenvalue,sigma] = fSimulate_AR2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_AR2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
else
    % MA(4)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_MA2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
end



%% set up the response process
[response, response_wo_noise, transfer_true, specCrossDensity_true] = fSimulate_response_ONB_v3(onb,response_process,response_process_shape, dataFullONB_padded,dataFull_padding,specDensity_true);
[response_test,~,~,~] = fSimulate_response_ONB_v3(onb,response_process,response_process_shape,dataFullONB_test_padded,dataFull_padding,specDensity_true);



%% visualization of {X_t}
figure('Position', [50,50, 1500, 800])
for t = 1:12
    subplot(3,4,t)
    plot(onb.gridSpace, dataFull(t,:))
    hold on
    plot( onb.gridSpace(censor.grid{t}) , censor.data{t},'x' )
    hold off
    
    dot_pr = mean( dataFull(t,:)' .* onb.onbMatrix * transfer_true.ONB(transfer_true.numOfLags + 1 ,:)' );
    
    integral = mean(dataFull(t,:));
    ylim([-4 4])
    xlim([0 1])
    title( ['t=',num2str(t)] )
end
suptitle('first 12 sparsely observed curves of \{X_t\}')
pause(0.1)


%% filters visualization
figure('Position', [100,100,1500, 800])
for lag = -3:3
    subplot(2,4,lag+4)
    plot( onb.gridSpace, onb.onbMatrix * transfer_true.ONB(transfer_true.numOfLags + lag + 1 ,:)' )
    ylim([-3 3])
    title("lag = "+num2str(lag) )
end
subplot(2,4,8)
plot( response(1:10), "-x" )
title("response t=1,..,10")
xlim([1 10])

suptitle("true filter coefficients; stretch of \{ Z_t \}")
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    3)      select the bandwidths for spatial smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% the mean is assumed to be zero
bw_mu = 1;
mu_est_smoother = zeros(onb.nGridSmoother,1);
mu_est_ONB = onb.smoothingMatrix_gridSmoother2ONB * mu_est_smoother;
mu_est_inSpace = onb.onbMatrix*mu_est_ONB;


%%    K-FOLD CV for the surface smoothing
[bw_r, bw_v, covLag0_smoother, covLag0_ridge_smoother, sigma2_est_loc_qv, sigma2_est_loc_lin] = fKCVforR_locQv(censor,onb,mu_est_inSpace);
covLag0_ONB = onb.smoothingMatrix_gridSmoother2ONB * covLag0_smoother * onb.smoothingMatrix_gridSmoother2ONB';
covLag0_inSpace = onb.onbMatrix * covLag0_ONB * onb.onbMatrix';

%% visualise
figure('Position', [150,150,1500, 800])
subplot(1,1,1)
surf( onb.gridSpace, onb.gridSpace, covLag0_inSpace )
title("estiamted lag-0 covariance kernel")


%%  K-FOLD CV for the cross-spectral density
[bw_c, crossCov0_smoother] = fKCVforCrossR(censor,onb,response,mu_est_inSpace);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    4)      estimate the spectral density and quantify the estimation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% estimate the spectral density
bw_F = bw_r; % the same smoothing parameter as for cov lag 0
specDensity_v7 = fEstimateSpecDensity_v7_2020(censor,onb,mu_est_inSpace,bw_F,nGridFreq);

% quantify the error of spec density estimation
v7_specDensity_error = fSpecDensityError( onb, specDensity_v7, specDensity_true );


%% visialise the spectral density
figure('Position', [200,200,1500, 800])
k_to_visualise = [500 600 700 800 900 1000];

real_minmax = [min(real(specDensity_v7.smoother(:))), max(real(specDensity_v7.smoother(:)))];
imag_minmax = [min(imag(specDensity_v7.smoother(:))), max(imag(specDensity_v7.smoother(:)))];


for ii = 1:length(k_to_visualise)
    k = k_to_visualise(ii);
    omega = specDensity_v7.gridFreq(k); % real omega: -pi to pi    
    
    m7 = onb.onbMatrix * squeeze(specDensity_v7.ONB_positified(k,:,:)) * onb.onbMatrix';    
    
    subplot(3,length(k_to_visualise),ii)
    surf( onb.gridSpace,onb.gridSpace, real(m7), 'EdgeColor','none')    
    title("omega = "+num2str(round(omega,2)))
    if ii == 1, zlabel("real part"); end
    zlim(real_minmax)
    
    subplot(3,length(k_to_visualise),ii+length(k_to_visualise))
    surf( onb.gridSpace,onb.gridSpace, imag(m7), 'EdgeColor','none')
    if ii == 1, zlabel("imag part"); end
    zlim(imag_minmax)
    
end

% draw norms over frequencies
subplot(3,1,3)
plot( specDensity_v7.gridFreq, specDensity_v7.trace_class_norms )
vline( specDensity_v7.gridFreq(k_to_visualise) )
ylabel("trace norm")
xlabel("frequency (highlighted the frequencies visualised above)")

suptitle("estimated spectral density \{F_\omega^X\}")





%% extract auto-covariances (integrate out the spectral density)
numOfLags_all = censor.nGridTime - 1;
v7_covLags = fEstimateCovLags(onb,specDensity_v7,numOfLags_all);

%% quantify the error of cov lag 0 estimation
covLag0_error = fCovLag0Error( onb, covsLags_true, covLag0_ONB );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    5)      functional recovery for the process X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if I want to get confidence bands
cofidence_bands = 0;

% use the sigma from the loc-qv estimation
s2 = sigma2_est_loc_qv;
if s2 < 1e-2, s2 = 1e-2; end
sigma_est = sqrt( s2 );

% kriging on the TRAINING sample
disp('functional data recovery')
[krigingX_dyna,krigingX_dyna_padded] = fKrigingX( censor, onb, v7_covLags, mu_est_ONB, sigma_est, 1, cofidence_bands, dataFull_padding );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    6)      estimate the CROSS-spectral density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('cross spectral density')
response_mean = 0;
specCrossDensity = fEstimateSpecCrossDensity_holdout_v5(censor,onb,response,response_mean,mu_est_inSpace,bw_c,nGridFreq,holdout); % bandwidth: bw_c



%% visualise the spec cross density
figure('Position', [250,250,1500, 800])

m = onb.onbMatrix * specCrossDensity.ONB';
subplot(2,2,1)
surf( specCrossDensity.gridFreq, onb.gridSpace, real(m) , 'EdgeColor','none')
xlabel("frequency")
ylabel("spatial dimension [0,1]")
zlabel("real part")
title("real part of F_\omega^{ZX}")


subplot(2,2,2)
surf( specCrossDensity.gridFreq, onb.gridSpace, imag(m) , 'EdgeColor','none')
xlabel("frequency")
ylabel("spatial dimension [0,1]")
zlabel("imag part")
title("imag part of F_\omega^{ZX}")

subplot(2,1,2)
norms = nan(specCrossDensity.nGridFreq,1);
for k = 1:specCrossDensity.nGridFreq
    norms(k) = norm(specCrossDensity.ONB(k,:)');
end

plot(specCrossDensity.gridFreq,norms)
xlabel("frequency")
ylabel("vector norm of F_\omega^{ZX}")
title("vector norm of F_\omega^{ZX} as a function of frequency")



suptitle("estimated spectral density \{F_\omega^{ZX}\}")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    7)      truncation regularisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% holdout cross-validation to determine the threshold
threshold_value_interval = [0.0001, 100];
[threshold_value] = fHoldoutCVforTransfer_trunc(threshold_value_interval, censor, onb, holdout, specDensity_v7, specCrossDensity, krigingX_dyna_padded, response, dataFull_padding);

%% fit the transfer functions
transfer_trunc = fEstimateTransfer_threshold_value(onb,specDensity_v7,specCrossDensity, threshold_value);        

% visualise
figure('Position', [300,300,1500, 800])

subplot(3,5,3*5)
plot( transfer_trunc.gridFreq, transfer_trunc.truncationEigen )%
title("truncation level")
xlabel("frequency")
ylabel("K_\omega")

for lag = -6:7
    subplot(3,5,lag+7)
    plot( onb.gridSpace, onb.onbMatrix * transfer_trunc.ONB(lag + transfer_trunc.numOfLags +1, :)' )
    if abs(lag) < transfer_true.numOfLags
        truth_ONB = transfer_true.ONB(lag + transfer_true.numOfLags +1, :)';
    else
        truth_ONB = zeros(onb.nBasis, 1);
    end
    hold on
    plot( onb.gridSpace, onb.onbMatrix * truth_ONB )
    hold off
    ylim([-3 3])
    title("lag="+num2str(lag))    
end

suptitle("Truncation regularisation; blue=estim, red=true")





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    8)      Tikhonov regularisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% holdout cross-validation to determine the threshold
tikhonovParameter_interval = [0.0001, 10];
[tikhonovParameter] = fHoldoutCVforTransfer_Tikh(tikhonovParameter_interval, censor, onb, holdout, specDensity_v7, specCrossDensity, krigingX_dyna_padded, response, dataFull_padding);

%% fit the transfer functions
transfer_Tikh = fEstimateTransfer_Tikhonov(onb,specDensity_v7,specCrossDensity, tikhonovParameter);


% visualise
figure('Position', [350,350,1500, 800])

for lag = -6:8
    subplot(3,5,lag+7)
    plot( onb.onbMatrix * transfer_Tikh.ONB(lag + transfer_Tikh.numOfLags +1, :)' )
    if abs(lag) < transfer_true.numOfLags
        truth_ONB = transfer_true.ONB(lag + transfer_true.numOfLags +1, :)';
    else
        truth_ONB = zeros(onb.nBasis, 1);
    end
    hold on
    plot( onb.onbMatrix * truth_ONB )
    hold off
    ylim([-3 3])
    title("lag="+num2str(lag))
end
suptitle("Tikhonov regularisation; blue=estim, red=true")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    9)      forecasting response process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% visualise
figure('Position', [400,400,1500, 800])

% estimation on the TRAINING sample
response_est_trunc = fKrigingZ(krigingX_dyna_padded,transfer_trunc,dataFull_padding);
response_est_Tikh = fKrigingZ(krigingX_dyna_padded,transfer_Tikh,dataFull_padding);

t_display = 1:100;

plot(response_est_trunc(t_display), "-x")
hold on
plot(response_est_Tikh(t_display), "-x")
plot(response(t_display), "-x")
hold off
legend("truncation","Tikhonov","truth")
title("snapshot of estimated response vs. ground truth")

