%{

This script loads the data saved variables in
"Wank_vis_estimated_dynamics_X.mat" and continues with the analysis of the
lagged regression dependence on the atmospheric electricity.

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}
clear all

addpath('../master')
addpath('../master/fdaM')

load('Wank_matlab_variables_after_a.mat')

%%
specDensity_v7_to_use = specDensity_v7s;
krigingX_dyna_to_use = krigingX_dyna_s;
krigingX_dyna_padded_to_use = krigingX_dyna_padded_s;

response_to_use = response;

%% visualise cov lag 0 kernel

subplot(1,1,1)
surf( covLag0_inSpace )
pause(.01)

%% diagonal
plot( diag(covLag0_inSpace) )
ylim([0,2000])


%% visialise the spectral density

for k = 1:100:specDensity_v7_to_use.nGridFreq
    omega = specDensity_v7_to_use.gridFreq(k); % real omega: -pi to pi        
    m7 = onb.onbMatrix * squeeze(specDensity_v7_to_use.ONB_positified(k,:,:)) * onb.onbMatrix';
    
    subplot(1,2,1)
    title('real')
    surf( real(m7), 'EdgeColor','none')        
    
    subplot(1,2,2)
    title('imag')
    surf( imag(m7), 'EdgeColor','none' )
    
    
    pause(.1)
    
end



%% some plots
% nDays = 3*5;
% offset = 80;
% for ii = 1:nDays
%     t = offset+ii;
%     subplot(3,5,ii)
%     
%     
%     confidence band
%     lower_dynamic     = onb.onbMatrix*krigingX_dyna.est(t,:)' - 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));
%     upper_dynamic     = onb.onbMatrix*krigingX_dyna.est(t,:)' + 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));       
%     z = fSCBDegras( squeeze(krigingX_dyna.var(t,:,:)), onb.onbMatrix );
%     lower_sim_dynamic = onb.onbMatrix*krigingX_dyna.est(t,:)' - z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));
%     upper_sim_dynamic = onb.onbMatrix*krigingX_dyna.est(t,:)' + z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna.var(t,:,:))*onb.onbMatrix')));
% 
%     fill([onb.gridSpace fliplr(onb.gridSpace)],[(lower_sim_dynamic') fliplr(upper_sim_dynamic')],'y', 'FaceAlpha', 0.4) % dynamic band (simultaneous)
% 
%     hold on
%     
%     observed data
%     plot( onb.gridSpace(censor.onb.grid24), Eat_daily(t,:) )    
%     nonzero_positions = find(Eat_indic(t,:));
%     scatter( onb.gridSpace(censor.onb.grid24(nonzero_positions)), Eat_daily(t, nonzero_positions), 'xb' )
% 
%     kriging visualisation
%     plot( onb.gridSpace, onb.onbMatrix * krigingX_dyna.est( t,: )' )
%     
%     hold off
%     title("day="+num2str(t))
%     
%     ylim([-100 300])
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    cross spectral density



%%  K-FOLD CV for the cross-spectral density
[bw_c, crossCov0_smoother] = fKCVforCrossR_NaNready(censor,onb,response_to_use,mu_est_inSpace);

%%
var_z = var( response_to_use, 'omitnan' );


crossCov0_inSpace = onb.onbMatrix * onb.smoothingMatrix_gridSmoother2ONB * crossCov0_smoother;

plot(crossCov0_inSpace / sqrt(var_z) ./ sqrt(diag(covLag0_inSpace))  )

%% mean of the response process
response_mean = mean(response_to_use, 'omitnan');


% the spec cross density
% holdout_partition = [0.7 0.85]; % 0 - 70% for fitting, 70% - 85% for CV, 85% - 100% for evaluation
holdout_end = 0.75;
% holdout_CV = 0.15;
% holdout_test = 0.15;

% bw_c_to_use = bw_c;
bw_c_to_use = 0.2;
E_specCrossDensity = fEstimateSpecCrossDensity_holdout_NaNready(censor,onb,response_to_use,...
    response_mean,mu_est_inSpace,bw_c_to_use,nGridFreq, [0 holdout_end]); % bandwidth: bw_c

% visualise the spectral cross density
m = onb.onbMatrix * E_specCrossDensity.ONB';
subplot(1,2,1)
surf( real(m), 'EdgeColor','none' )
subplot(1,2,2)
surf( imag(m), 'EdgeColor','none' )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% truncation and Tikhonov regularization to fit the filter coefficients


% holdout cross-validation to determine the threshold
threshold_value_interval = [0.0001, 100];
[threshold_value] = fHoldoutCVforTransfer_w_mean_trunc(threshold_value_interval, censor, onb, [holdout_end 1], specDensity_v7, E_specCrossDensity, krigingX_dyna_padded_to_use, response_to_use, dataFull_padding,mu_est_ONB,response_mean);
tikhonovParameter_interval = [10 4000];
[tikhonovParameter] = fHoldoutCVforTransfer_w_mean_Tikh(tikhonovParameter_interval, censor, onb, [holdout_end 1], specDensity_v7, E_specCrossDensity, krigingX_dyna_padded_to_use, response_to_use, dataFull_padding,mu_est_ONB,response_mean);

threshold_value
tikhonovParameter

%% fit the transfer functions
threshold_value = 10;
tikhonovParameter = 35;
trunc_E_transfer = fEstimateTransfer_threshold_value(onb,specDensity_v7_to_use,E_specCrossDensity, threshold_value);        
Tikh_E_transfer = fEstimateTransfer_Tikhonov(onb,specDensity_v7_to_use,E_specCrossDensity, tikhonovParameter);        

% visualise
% suptitle("Truncation")
% subplot(3,5,20)
% plot( transfer.truncationEigen )%
% title("truncation level")
for lag = -7:7    
    subplot(3,5,-lag+8) % 4*5
    plot( onb.gridSpace, onb.onbMatrix * trunc_E_transfer.ONB(lag + trunc_E_transfer.numOfLags +1, :)' )    
    hold on
    plot( onb.gridSpace, onb.onbMatrix * Tikh_E_transfer.ONB(lag + Tikh_E_transfer.numOfLags +1, :)' )    
    hold off
    ylim([-1.5 1.5])
%     ylim([-1 1]*0.05)
    title("lag="+num2str(lag))
    legend({"trunc","Tikh"})
end
pause(.1)


%% plot with missing data
plot(response_to_use)
vline( find(~response_available) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualisation of prediction

trunc_E_response_est = fKrigingZ_nonzeromean(krigingX_dyna_padded_to_use,trunc_E_transfer, mu_est_ONB, response_mean ,dataFull_padding);
Tikh_E_response_est = fKrigingZ_nonzeromean(krigingX_dyna_padded_to_use,Tikh_E_transfer, mu_est_ONB, response_mean ,dataFull_padding);

% where to start
test_start = ceil( censor.nGridTime * holdout_end + 0.001);

t_view = test_start:censor.nGridTime;

subplot(2,1,1)
plot( response_to_use(t_view) )
hold on
plot( trunc_E_response_est(t_view) )
plot( Tikh_E_response_est(t_view) )
hline(response_mean)
hold off
legend(["truth","est.trunc","est.Tikh"])
title("fit")
% ylim([-20 100])

disp( "mean of response  "+num2str(mean(response_to_use, 'omitnan')) )
disp( "mean of est.trunc "+num2str(mean(trunc_E_response_est, 'omitnan')) )
disp( "mean of est.Tikh  "+num2str(mean(Tikh_E_response_est, 'omitnan')) )

x_vline = 75;
vline(x_vline)
% hline(88.0629)

% plot residuals
subplot(2,1,2)
plot( response_to_use(t_view)- response_mean*ones(size(Tikh_E_response_est(t_view)))'  )
hold on
plot( response_to_use(t_view)-Tikh_E_response_est(t_view)'  )
plot( response_to_use(t_view)-trunc_E_response_est(t_view)'  )
hline(0)
hold off
legend(["mean.pred","est.trunc","est.Tikh"])
title("residuals")
% ylim([-55 65])
vline(x_vline)

t = t_view(1)-1+x_vline;
response_to_use(t)
trunc_E_response_est(t)
Tikh_E_response_est(t)


%% MSE, R2 train partition
% 
% mse_test_trunc = fKrigingZ_holdout_MSE_NaNready( response_est, response', [0 holdout_partition(1)])
% mse_test_Tikh  = fKrigingZ_holdout_MSE_NaNready( response_est_Tikh, response', [0 holdout_partition(1)])
% 
% ss_tot = fKrigingZ_holdout_MSE_NaNready( response_mean*ones(size(response_est)), response', [0 holdout_partition(1)])
% 
% r2_trunc = 1 - mse_test_trunc/ss_tot
% r2_Tikh = 1 - mse_test_Tikh/ss_tot

%% MSE on the test set

mse_test_trunc = fKrigingZ_holdout_MSE_NaNready( trunc_E_response_est, response_to_use', [holdout_end 1])
mse_test_Tikh  = fKrigingZ_holdout_MSE_NaNready( Tikh_E_response_est, response_to_use', [holdout_end 1])

ss_tot = fKrigingZ_holdout_MSE_NaNready( response_mean*ones(size(trunc_E_response_est)), response_to_use', [holdout_end 1])

r2_trunc = 1 - mse_test_trunc/ss_tot
r2_Tikh = 1 - mse_test_Tikh/ss_tot



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%      mean of Kriged process
% 
% x = mean( onb.onbMatrix * krigingX_dyna_padded_to_use.est'  , 2);
% size(x)
% plot(x)
% hold on
% plot(mu_est_inSpace)
% hold off
% legend('kriged mean','estimated daily mean')
% 
% %% calculate theoretical response mean after applying the response coefficients
% x = 0;
% for lag = -transfer.numOfLags:transfer.numOfLags
%     x = x + transfer.ONB(lag + transfer.numOfLags +1, :) * mu_est_ONB;
% end
% x



%% with or withouth mean (substracted mu_est_ONB)

% with mean
% adjustment = zeros(size(mu_est_ONB));
% Eat_daily_centred = Eat_daily;

% without mean
adjustment = mu_est_ONB;
Eat_daily_centred = Eat_daily - mean(Eat_daily,1,'omitnan');



nDays = 3*5;
offset = 80;
day_center = 88;
contribution_total = 0;
contributions = [];
for ii = 1:(nDays)
    t = offset+ii;
    
    if (t >= 1) && (t <= censor.nGridTime)
        subplot(3,5,ii)    

        % confidence band
        lower_dynamic     = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) - 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));
        upper_dynamic     = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) + 1.96*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));       
        z = fSCBDegras( squeeze(krigingX_dyna_to_use.var(t,:,:)), onb.onbMatrix );
        lower_sim_dynamic = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) - z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));
        upper_sim_dynamic = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) + z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));

        fill([onb.gridSpace fliplr(onb.gridSpace)],[(lower_sim_dynamic') fliplr(upper_sim_dynamic')],'y', 'FaceAlpha', 0.4) % dynamic band (simultaneous)

        hold on

        % observed data
%         plot( onb.gridSpace(censor.onb.grid24), Eat_daily_centred(t,:) )    
        nonzero_positions = find(Eat_indic(t,:));
        scatter( onb.gridSpace(censor.onb.grid24(nonzero_positions)), Eat_daily_centred(t, nonzero_positions), 'xb' )

        % kriging visualisation
        plot( onb.gridSpace, onb.onbMatrix * (krigingX_dyna.est( t,: )'-adjustment) )
        
        ylim([-150 200])

        hold off
    end
    
    % contribution from day "day_center"
%     lag_real = -day_center + t;
    lag_real = +day_center - t;
    lag_index=lag_real + trunc_E_transfer.numOfLags+1;
    if abs(lag_real)<= trunc_E_transfer.numOfLags
        contribution = trunc_E_transfer.ONB(lag_index,:) * ( krigingX_dyna_padded_to_use.est( t+dataFull_padding ,:)' - mu_est_ONB )    ;
        contributions = [contributions contribution];
        contribution_total = contribution_total + contribution;
        disp("t="+num2str(t)+", lag_real="+num2str(lag_real)+", lag_index="+num2str(lag_index)+", contribution="+num2str(contribution))
    
    end
    title("lag="+num2str(lag_real)+", day="+num2str(t)+", contr="+num2str(contribution))
end

suptitle("mean centred predictions, contributions for day 88")

contributions = flipud(contributions');

%
disp("***************************************************")
disp("result from visualisation = "+num2str(contribution_total + response_mean))
disp("result from estimfunction = "+num2str(trunc_E_response_est(day_center)))


