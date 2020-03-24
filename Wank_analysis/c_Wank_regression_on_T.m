%{

This script continues with the analysis of the
lagged regression dependence on the temperature.

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}



%% fully function TS for temperature

% extract the days
T_daily = nan(nGridTime, 24);

day = 0;
DOY_last = -9999;
for t = 1:height(dataFull)
    
    DOY = dataFull.DOY(t)+1;
    if DOY ~= DOY_last, day = day+1; end
    DOY_last = DOY;
    hour = dataFull.hour(t)+1;
    
    T_daily(day,hour) = dataFull.temp_no_seasonality(t);
end

%% visualisation

nDays = 3*5;
offset = 80;
for ii = 1:nDays
    t = offset+ii;
    subplot(3,5,ii)
    
    
    % confidence band
    lower_dynamic     = onb.onbMatrix*krigingX_dyna.est(t,:)' - 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));
    upper_dynamic     = onb.onbMatrix*krigingX_dyna.est(t,:)' + 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));       
    z = fSCBDegras( squeeze(krigingX_dyna.var(t,:,:)), onb.onbMatrix );
    lower_sim_dynamic = onb.onbMatrix*krigingX_dyna.est(t,:)' - z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));
    upper_sim_dynamic = onb.onbMatrix*krigingX_dyna.est(t,:)' + z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));

    fill([onb.gridSpace fliplr(onb.gridSpace)],[(lower_sim_dynamic') fliplr(upper_sim_dynamic')],'y', 'FaceAlpha', 0.4) % dynamic band (simultaneous)

    hold on
    
    % observed data
%     plot( onb.gridSpace(censor.onb.grid24), Eat_daily(t,:) )    
    nonzero_positions = find(Eat_indic(t,:));
    scatter( onb.gridSpace(censor.onb.grid24(nonzero_positions)), Eat_daily(t, nonzero_positions), 'xb' )

    % kriging visualisation
    plot( onb.gridSpace, onb.onbMatrix * krigingX_dyna.est( t,: )' )
    
    % temperature
    plot( onb.gridSpace(censor.onb.grid24), T_daily(t,:) * 10)
    
    hold off
    title("day="+num2str(t)+"; vis="+num2str(round(response(t))))
    
    ylim([-100 300])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              fully functional analysis of T_daily
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% holdout_partition = [0.70, 0.85]; % 0 - 70% for fitting, 70% - 85% for CV, 85% - 100% for evaluation

%% express the data in ONB
dataFullONB = nan(nGridTime, onb.nBasis);
T_is_available = ones(nGridTime);

dataFullONB(1,:) =  onb.smoothingMatrix_gridSmoother2ONB * [T_daily(1,1), T_daily(1,:) ]';
if sum(isnan(dataFullONB(1,:))) > 0, T_is_available(1) = 0; end
for t = 2:nGridTime    
    dataFullONB(t,:) =  onb.smoothingMatrix_gridSmoother2ONB * [T_daily(t-1,24), T_daily(t,:) ]';
    if sum(isnan(dataFullONB(t,:))) > 0, T_is_available(t) = 0; end
end

%
nDays = 2*3;
offset = 0;
for ii = 1:nDays
    t = offset+ii;
    subplot(2,3,ii)
    
    plot( onb.gridSpace, onb.onbMatrix * dataFullONB(t,:)')
    hold on
    plot( onb.gridSpace(censor.onb.grid24+5), T_daily(t,:), '-x')
    hold off
end

%% missingness plot
plot( mean(T_daily,2) )
vline(find(~T_is_available))


%% mean function
T_mean = mean(dataFullONB,'omitnan')';
plot( onb.onbMatrix * T_mean )


%% calculate response mean
response_mean = mean(response_to_use,'omitnan');


%% prepare the estimates of the autocovariance operators
holdout_start = ceil( nGridTime * holdout_end + 0.001);

num_of_lags = nGridTime;
cov_estim = nan(num_of_lags,onb.nBasis,onb.nBasis);
crosscov_estim = nan(num_of_lags,onb.nBasis);
crosscov_estim_minus_lags = nan(num_of_lags,onb.nBasis);
for lag_real = 0:(num_of_lags-1)
    lag_index = lag_real + 1;
    
    % covariances
    n_contributions = 0;
    sums = zeros(onb.nBasis);
    for t2 = 1:(nGridTime-lag_real)
        t1 = t2 + lag_real;
        
        m = (dataFullONB( t1, :)-T_mean')' * (dataFullONB( t2, :)-T_mean');
        if sum(isnan(m)) == 0
            n_contributions = n_contributions + 1;
            sums = sums + m;
        end
    end
    cov_estim(lag_index,:,:) = sums / n_contributions;
    
    % cross-covariances    
    n_contributions = 0;
    sums = zeros(1,onb.nBasis);
    for t2 = 1:(holdout_start-1-lag_real)
        t1 = t2 + lag_real;
        
        m = (response_to_use( t1 )-response_mean) * (dataFullONB( t2 ,:)-T_mean');
        if sum(isnan(m)) == 0
            n_contributions = n_contributions + 1;
            sums = sums + m;
        end
    end
    crosscov_estim(lag_index,:) = sums / n_contributions;
    
    % cross-covariances, minus lags
    n_contributions = 0;
    sums = zeros(1,onb.nBasis);
    for t2 = (lag_real+1):(holdout_start-1)
        t1 = t2 - lag_real;
        
        m = (response_to_use( t1 )-response_mean) * (dataFullONB( t2 ,:)-T_mean');
        if sum(isnan(m)) == 0
            n_contributions = n_contributions + 1;
            sums = sums + m;
        end
    end
    crosscov_estim_minus_lags(lag_index,:) = sums / n_contributions;
end

%% visualise - covariance surface
subplot(1,1,1)
lag_real = 0;
surf(onb.onbMatrix * squeeze(cov_estim(lag_real + 1,:,:)) * onb.onbMatrix');

%% diag
variances = diag(onb.onbMatrix * squeeze(cov_estim(lag_real + 1,:,:)) * onb.onbMatrix');
plot(variances)

%% correlation surface
lag_real = 0;
surf(onb.onbMatrix * squeeze(cov_estim(lag_real + 1,:,:)) * onb.onbMatrix' ./ sqrt(variances * variances') );
zlim([0 1])

%% cross covariance
lag_real = 0;
response_var = var(response_to_use,'omitnan');
plot( onb.onbMatrix * crosscov_estim(lag_real+1,:)'  ./ sqrt(variances) / sqrt(response_var) )
ylim([0 1])

%% estimate the spectral density and specCrossDensity
nGridFreq = 1000;
T_specDensity = [];
T_specDensity.nGridFreq = nGridFreq;
T_specDensity.gridFreq = linspace(-pi,pi,T_specDensity.nGridFreq);
q_bartlett = floor(sqrt(nGridTime));
T_specDensity.q_bartlett = q_bartlett;
T_specDensity.ONB = nan(nGridFreq, onb.nBasis, onb.nBasis);

T_specCrossDensity = [];
T_specCrossDensity.nGridFreq = nGridFreq;
T_specCrossDensity.gridFreq = linspace(-pi,pi, T_specCrossDensity.nGridFreq);
T_specCrossDensity.ONB = nan(T_specCrossDensity.nGridFreq, onb.nBasis);
T_specCrossDensity.q_bartlett = q_bartlett;

for k = 1:nGridFreq    
    omega = T_specDensity.gridFreq(k);
    % zero lag
    T_specDensity.ONB(k,:,:) = 1/(2*pi)*squeeze(cov_estim(1,:,:));
    T_specCrossDensity.ONB(k,:) = 1/(2*pi)*crosscov_estim(1,:);
    % non-zero lags
    for lag_real = 1:(q_bartlett-1)
        lag_index = lag_real + 1;
        T_specDensity.ONB(k,:,:) = squeeze(T_specDensity.ONB(k,:,:)) +... % positive lag
            1/(2*pi)*(1 - lag_real/q_bartlett)*squeeze(cov_estim(lag_index,:,:)) * exp(-1i*omega*lag_real);
        T_specDensity.ONB(k,:,:) = squeeze(T_specDensity.ONB(k,:,:)) +... % negative lag => transpose the autocovariance
            1/(2*pi)*(1 - lag_real/q_bartlett)*squeeze(cov_estim(lag_index,:,:))' * exp(+1i*omega*lag_real);
        T_specCrossDensity.ONB(k,:) = T_specCrossDensity.ONB(k,:) +...
            1/(2*pi)*(1 - lag_real/q_bartlett)*crosscov_estim(lag_index,:) * exp(-1i*omega*lag_real);
        T_specCrossDensity.ONB(k,:) = T_specCrossDensity.ONB(k,:) +...
            1/(2*pi)*(1 - lag_real/q_bartlett)*crosscov_estim_minus_lags(lag_index,:) * exp(+1i*omega*lag_real);
    end
    
end

% positify the spectral densities, truncate the space corresponding to negative eigenvalues
T_specDensity.ONB_positified = zeros(nGridFreq, onb.nBasis, onb.nBasis);
for k = 1:nGridFreq
    m = squeeze(T_specDensity.ONB(k,:,:));
    m = (m+m')/2; % make sure it's self-adjoint
    [V,D] = eig(m); % i.e. the decomposition:  m = V*D*V'
    diagD = real(diag(D)); % note that the harmonic eigenvalues must be always real
    diagD( diagD< 0 ) = 0; % truncate those elements that are less than (numerical) zero
    m = V*diag(diagD)*V';
    T_specDensity.ONB_positified(k,:,:) = (m+m')/2; % again, make sure it's self-adjoint  
end


%% traces
traces = nan(T_specDensity.nGridFreq,1);

for k = 1:T_specDensity.nGridFreq
    m7 = onb.onbMatrix * squeeze(T_specDensity.ONB_positified(k,:,:)) * onb.onbMatrix';
    traces(k) = trace(m7);
end

subplot(1,1,1)
plot(real(traces))

%% visualise specDensity

for k = 1:100:T_specDensity.nGridFreq
    omega = T_specDensity.gridFreq(k); % real omega: -pi to pi        
    m7 = onb.onbMatrix * squeeze(T_specDensity.ONB_positified(k,:,:)) * onb.onbMatrix';
    
    subplot(1,2,1)
    title('real')
    surf( real(m7), 'EdgeColor','none')        
    
    subplot(1,2,2)
    title('imag')
    surf( imag(m7), 'EdgeColor','none' )
    
    pause(.1)
    
end

%% visualise specCrossDensity

subplot(1,2,1)
surf(real( onb.onbMatrix * T_specCrossDensity.ONB' ), 'edgecolor','none')
hold on
subplot(1,2,2)
surf(imag( onb.onbMatrix * T_specCrossDensity.ONB' ), 'edgecolor','none')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regularization

% here I'm padding with mean function
dataFullONB_padded_w_zeros=[];
dataFullONB_padded_w_zeros.est = [nan(dataFull_padding,onb.nBasis); dataFullONB; nan(dataFull_padding,onb.nBasis)];

% replace NaN with the mean function
for t = 1:(nGridTime + 2*dataFull_padding)
    if sum(isnan(dataFullONB_padded_w_zeros.est(t,:))) > 0
        dataFullONB_padded_w_zeros.est(t,:) = T_mean' ;
    end
end


% CV for regularization

% % holdout cross-validation to determine the threshold
% threshold_value_interval = [0.01, 1];
% [threshold_value] = fHoldoutCVforTransfer_trunc(threshold_value_interval, censor, onb, holdout_partition, T_specDensity, T_specCrossDensity, dataFullONB_padded_w_zeros, response_to_use, dataFull_padding,T_mean,response_mean);
% tikhonovParameter_interval = [0.01 1];
% [tikhonovParameter] = fHoldoutCVforTransfer_Tikh(tikhonovParameter_interval, censor, onb, holdout_partition, T_specDensity,T_specCrossDensity, dataFullONB_padded_w_zeros, response_to_use, dataFull_padding,T_mean,response_mean);
% 
% threshold_value
% tikhonovParameter

%% fit and visualise
threshold_value = 0.2
tikhonovParameter = 0.8
transfer = fEstimateTransfer_threshold_value(onb,T_specDensity,T_specCrossDensity, threshold_value);        
transfer_Tikh = fEstimateTransfer_Tikhonov(onb,T_specDensity,T_specCrossDensity, tikhonovParameter);        


for lag = -7:7    
    subplot(3,5,-lag+8) % 4*5
    plot( onb.gridSpace, onb.onbMatrix * transfer.ONB(lag + transfer.numOfLags +1, :)' )    
    hold on
    plot( onb.gridSpace, onb.onbMatrix * transfer_Tikh.ONB(lag + transfer_Tikh.numOfLags +1, :)' )    
    hold off
%     ylim([-1.5 1.5])
    ylim([-1 1]*5)
    title("lag="+num2str(lag))
    legend({"trunc","Tikh"})
end
pause(.1)

% truncation levels
subplot(3,5,15)
plot( transfer.truncationEigen )%
title("truncation level")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualisation of prediction

trunc_T_response_est = fKrigingZ_nonzeromean(dataFullONB_padded_w_zeros,transfer, T_mean, response_mean ,dataFull_padding);
Tikh_T_response_est = fKrigingZ_nonzeromean(dataFullONB_padded_w_zeros,transfer_Tikh, T_mean, response_mean ,dataFull_padding);

% where to start
test_start = ceil( censor.nGridTime * holdout_end + 0.001);

t_view = test_start:censor.nGridTime;

subplot(2,1,1)
plot( response_to_use(t_view) )
hold on
plot( trunc_T_response_est(t_view) )
plot( Tikh_T_response_est(t_view) )
hline(response_mean)
hold off
legend(["truth","est.trunc","est.Tikh"])
title("fit")
% ylim([-20 100])

disp( "mean of response  "+num2str(mean(response_to_use, 'omitnan')) )
disp( "mean of est.trunc "+num2str(mean(trunc_T_response_est, 'omitnan')) )
disp( "mean of est.Tikh  "+num2str(mean(Tikh_T_response_est, 'omitnan')) )

x_vline = 98;
vline(x_vline)
% hline(88.0629)

% plot residuals
subplot(2,1,2)
plot( response_to_use(t_view)- response_mean*ones(size(Tikh_T_response_est(t_view)))'  )
hold on
plot( response_to_use(t_view)-Tikh_T_response_est(t_view)'  )
plot( response_to_use(t_view)-trunc_T_response_est(t_view)'  )
hline(0)
hold off
legend(["mean.pred","est.trunc","est.Tikh"])
title("residuals")
% ylim([-55 65])
vline(x_vline)
% 
t = t_view(1)-1+x_vline;
% response_to_use(t)
% response_est(t)
% response_est_Tikh(t)


% MSE on the test set

mse_test_trunc = fKrigingZ_holdout_MSE_NaNready( trunc_T_response_est, response_to_use', [holdout_end 1])
mse_test_Tikh  = fKrigingZ_holdout_MSE_NaNready( Tikh_T_response_est, response_to_use', [holdout_end 1])

ss_tot = fKrigingZ_holdout_MSE_NaNready( response_mean*ones(size(trunc_T_response_est)), response_to_use', [holdout_end 1])

r2_trunc = 1 - mse_test_trunc/ss_tot
r2_Tikh = 1 - mse_test_Tikh/ss_tot
