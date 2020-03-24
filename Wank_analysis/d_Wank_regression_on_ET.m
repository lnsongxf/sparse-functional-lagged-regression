%{

This script continues with the analysis of the
joint lagged regression dependence on the atmoshperic electricity and temperature.
The files "b_Wank_regression_on_E.m" and "c_Wank_regression_on_E.m" should have been run before.

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}


save_plots = 1;

%%

ET_specCrossDensity = [];
ET_specCrossDensity.nGridFreq = nGridFreq;
ET_specCrossDensity.gridFreq = linspace(-pi,pi, ET_specCrossDensity.nGridFreq);
ET_specCrossDensity.smoother = nan(ET_specCrossDensity.nGridFreq, onb.nGridSmoother, onb.nGridSmoother);
ET_specCrossDensity.ONB = nan(ET_specCrossDensity.nGridFreq, onb.nBasis, onb.nBasis);
ET_specCrossDensity.q_bartlett = q_bartlett;

for ii = 1:onb.nGridSmoother
    x_smoother = onb.gridSmoother(ii);
    
    response_on_T = (onb.onbMatrix(onb.gridSpace2Smoother_indx(ii),:) * dataFullONB');
    response_on_T_mean = mean(response_on_T, 'omitnan');
    
    specCrossDensity = fEstimateSpecCrossDensity_holdout_NaNready(censor,onb,response_on_T,response_on_T_mean,mu_est_inSpace,bw_c_to_use,nGridFreq, [0 holdout_end]);
    
    ET_specCrossDensity.smoother(:,ii,:) = specCrossDensity.smoother;
end

%% convert to ONB
for k = 1:ET_specCrossDensity.nGridFreq
    ET_specCrossDensity.ONB(k,:,:) = onb.smoothingMatrix_gridSmoother2ONB * squeeze(ET_specCrossDensity.smoother(k,:,:)) * onb.smoothingMatrix_gridSmoother2ONB';
end


%% visialise the T_E_specCrossDensity

for k = 1:100:ET_specCrossDensity.nGridFreq
    omega = ET_specCrossDensity.gridFreq(k); % real omega: -pi to pi
    m = onb.onbMatrix * squeeze(ET_specCrossDensity.ONB(k,:,:)) * onb.onbMatrix';
    
    subplot(1,2,1)
    title('real')
    surf( real(m), 'EdgeColor','none')
    
    subplot(1,2,2)
    title('imag')
    surf( imag(m), 'EdgeColor','none' )
    
    
    pause(.01)
    
end

%% visualise the join spectral density


for k = 1:100:ET_specCrossDensity.nGridFreq
    omega = ET_specCrossDensity.gridFreq(k); % real omega: -pi to pi
    m_E = onb.onbMatrix * squeeze(specDensity_v7_to_use.ONB(k,:,:)) * onb.onbMatrix';
    m_T = onb.onbMatrix * squeeze(T_specDensity.ONB(k,:,:)) * onb.onbMatrix' * 100;
    m_ET = onb.onbMatrix * squeeze(ET_specCrossDensity.ONB(k,:,:)) * onb.onbMatrix' * 10;
    
    subplot(1,2,1)
    surf( onb.gridSpace,onb.gridSpace, real(m_E), 'EdgeColor','none')
    hold on
    surf( onb.gridSpace+1,onb.gridSpace+1, real(m_T), 'EdgeColor','none')
    surf( onb.gridSpace,onb.gridSpace+1, real(m_ET), 'EdgeColor','none')
    surf( onb.gridSpace+1,onb.gridSpace, real(m_ET'), 'EdgeColor','none')
    hold off
    title("real, omega="+num2str(omega))
    zlim([-1 1]*2000)
    
    subplot(1,2,2)
    surf( onb.gridSpace,onb.gridSpace, imag(m_E), 'EdgeColor','none')
    hold on
    surf( onb.gridSpace+1,onb.gridSpace+1, imag(m_T), 'EdgeColor','none')
    surf( onb.gridSpace,onb.gridSpace+1, imag(m_ET), 'EdgeColor','none')
    surf( onb.gridSpace+1,onb.gridSpace, imag(m_ET'), 'EdgeColor','none')
    hold off
    title('imag')
    zlim([-1 1]*300)
    
    pause(.001)
    
end

% %% correlations of E and T
% disp("calculating E_covLags")
% E_covLags = fEstimateCovLags(onb,specDensity_v7s,numOfLags_all);
% disp("calculating T_covLags")
% T_covLags = fEstimateCovLags(onb,T_specDensity,numOfLags_all);
%
% %% visualise
% subplot(2,2,1)
% surf( squeeze(E_covLags.corr_inSpace(1,:,:)), 'EdgeColor', 'none' )
% title("E corr 0")
%
% subplot(2,2,2)
% surf( squeeze(E_covLags.corr_inSpace(2,:,:)), 'EdgeColor', 'none' )
% title("E corr 1")
%
% subplot(2,2,3)
% surf( squeeze(T_covLags.corr_inSpace(1,:,:)), 'EdgeColor', 'none' )
% title("T corr 0")
%
% subplot(2,2,4)
% surf( squeeze(T_covLags.corr_inSpace(2,:,:)), 'EdgeColor', 'none' )
% title("T corr 1")
%
%
%
% %% cross correlation analysis
%
% E_cov_diag = diag(onb.onbMatrix* squeeze(E_covLags.ONB(1,:,:)) * onb.onbMatrix');
% T_cov_diag = diag(onb.onbMatrix* squeeze(E_covLags.ONB(1,:,:)) * onb.onbMatrix');
%
% ET_crossCovLags = [];
% % integrate spectral density to get covariances
% ET_crossCovLags.corr_numOfLags = 10;
% ET_crossCovLags.cov_inSpace = zeros(2*ET_crossCovLags.corr_numOfLags+1,onb.nGridSpace,onb.nGridSpace);
% ET_crossCovLags.corr_inSpace = zeros(2*ET_crossCovLags.corr_numOfLags+1,onb.nGridSpace,onb.nGridSpace);
%
% for lag_real=-ET_crossCovLags.corr_numOfLags:ET_crossCovLags.corr_numOfLags
%     lag_index = lag_real + ET_crossCovLags.corr_numOfLags + 1;
%
%     m = zeros(onb.nBasis);
%     for k=1:ET_specCrossDensity.nGridFreq % Fourier the positified spectral density
%         omega = ET_specCrossDensity.gridFreq(k);
%         m = m + squeeze(ET_specCrossDensity.ONB(k,:,:)) * exp(1i*omega*lag_real) / ET_specCrossDensity.nGridFreq*2*pi;
%     end
%
%     m = onb.onbMatrix * real(m) * onb.onbMatrix';
%
%     ET_crossCovLags.cov_inSpace(lag_index,:,:) = m;
%     ET_crossCovLags.corr_inSpace(lag_index,:,:) = m ./ sqrt( E_cov_diag*T_cov_diag');
% end
%
% %%
% for lag_real = -3:3
%     lag_index = lag_real + ET_crossCovLags.corr_numOfLags + 1;
%     subplot(2,4,lag_real+4)
%     surf( squeeze(ET_crossCovLags.corr_inSpace(lag_index,:,:)), 'EdgeColor','none' )
%     title("lag = "+num2str(lag_real))
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tikhonov regularisation

tikhonovParameter_1 = 35;
tikhonovParameter_2 = 0.8;

[Tikh_transfer_1,Tikh_transfer_2] = fEstimateTransfer_Tikhonov_dual(onb,...
    specDensity_v7_to_use,T_specDensity,ET_specCrossDensity,E_specCrossDensity,T_specCrossDensity,...
    tikhonovParameter_1, tikhonovParameter_2);

truncationParameter_1 = 10;
truncationParameter_2 = 0.3;

[trunc_transfer_1,trunc_transfer_2] = fEstimateTransfer_truncation_dual(onb,...
    specDensity_v7_to_use,T_specDensity,ET_specCrossDensity,E_specCrossDensity,T_specCrossDensity,...
    truncationParameter_1, truncationParameter_2);

% visualise

for lag_real = -7:7
    subplot(3,5,-lag_real+8)
    plot( onb.gridSpace, onb.onbMatrix * trunc_transfer_1.ONB(lag_real + trunc_transfer_1.numOfLags +1, :)' )
    hold on
    plot( onb.gridSpace, onb.onbMatrix * trunc_transfer_2.ONB(lag_real + trunc_transfer_2.numOfLags +1, :)' )
    plot( onb.gridSpace, onb.onbMatrix * Tikh_transfer_1.ONB(lag_real + Tikh_transfer_1.numOfLags +1, :)' )
    plot( onb.gridSpace, onb.onbMatrix * Tikh_transfer_2.ONB(lag_real + Tikh_transfer_2.numOfLags +1, :)' )
    hold off
    hline(0)
    %     ylim([-1.5 1.5])
    ylim([-1 1]*5)
    title("lag="+num2str(lag_real))
    legend({"trunc.E","trunc.T","Tikh.E","Tikh.T"})
end
pause(.1)

%% visualise only filters E and T on different rows


if save_plots, h_fig = figure('Position', [200,200,800, 400]); end


for lag_real = -3:3
    subplot(2,7,-lag_real+4)
    plot( 24*onb.gridSpace, onb.onbMatrix * trunc_transfer_1.ONB(lag_real + trunc_transfer_1.numOfLags +1, :)', "--k" )
    hold on
    plot( 24*onb.gridSpace, onb.onbMatrix * Tikh_transfer_1.ONB(lag_real + Tikh_transfer_1.numOfLags +1, :)', "-k" )
    hold off
    hline(0, ":k")
    ylim([-1.5 1.5])
    title("B^{(E)}_{"+num2str(lag_real)+"}" )
    xticks([0 8 16 24])
    
    if lag_real == 3
        ylabel({"estimated filter coefficients";"for electricity TS"})
    end
    
    subplot(2,7,-lag_real+4+7)
    plot( 24*onb.gridSpace, onb.onbMatrix * trunc_transfer_2.ONB(lag_real + trunc_transfer_2.numOfLags +1, :)', "--k"  )
    hold on
    plot( 24*onb.gridSpace, onb.onbMatrix * Tikh_transfer_2.ONB(lag_real + Tikh_transfer_2.numOfLags +1, :)', "-k" )
    hold off
    hline(0, ":k")
    ylim([-1 1]*6)
    title("B^{(\tau)}_{"+num2str(lag_real)+"}" )
    xticks([0 8 16 24])
    
    if lag_real == 3
        ylabel({"estimated filter coefficients";"for temperature TS"})
    end
end
pause(.1)

if save_plots
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, 'figures_paper/Wank_filters.pdf','-dpdf','-r0')
    close(h_fig)
end


%% prediction

Tikh_ET_response_est = fKrigingZ_dual( krigingX_dyna_padded_to_use, dataFullONB_padded_w_zeros,...
    Tikh_transfer_1,Tikh_transfer_2, mu_est_ONB, T_mean, response_mean ,dataFull_padding);
trunc_ET_response_est = fKrigingZ_dual( krigingX_dyna_padded_to_use, dataFullONB_padded_w_zeros,...
    trunc_transfer_1,trunc_transfer_2, mu_est_ONB, T_mean, response_mean ,dataFull_padding);

% where to start
test_start = ceil( censor.nGridTime * holdout_end + 0.001);

t_view = test_start:censor.nGridTime;
% t_view = 847:849

subplot(1,1,1)
plot( response_to_use(t_view) )
hold on
plot( Tikh_ET_response_est(t_view) )
plot( Tikh_E_response_est(t_view) )
plot( Tikh_T_response_est(t_view) )
hline(response_mean)
hold off
legend(["truth","ET.Tikh","E.Tikh","T.Tikh"])
title("fit")
ylim([-20 100])
% xlim([0 50])

% % plot residuals
% subplot(2,1,2)
% plot( response_to_use(t_view)- response_mean*ones(size(ET_response_est(t_view)))'  )
% hold on
% plot( response_to_use(t_view)-ET_response_est(t_view)'  )
% plot( response_to_use(t_view)-E_response_est_Tikh(t_view)'  )
% plot( response_to_use(t_view)-T_response_est_Tikh(t_view)'  )
% hline(0)
% hold off
% legend(["mean.pred","ET.Tikh","E.Tikh","T.Tikh"])
% title("residuals")



% MSE on the test set

% % log response
% ET_mse_test_Tikh  = fKrigingZ_holdout_MSE_NaNready( exp(ET_response_est), exp(response_to_use)', [holdout_end 1]);
% E_mse_test_Tikh   = fKrigingZ_holdout_MSE_NaNready( exp(E_response_est_Tikh), exp(response_to_use)', [holdout_end 1]);
% T_mse_test_Tikh   = fKrigingZ_holdout_MSE_NaNready( exp(T_response_est_Tikh), exp(response_to_use)', [holdout_end 1]);
% ss_tot = fKrigingZ_holdout_MSE_NaNready( exp(response_mean)*ones(size(ET_response_est)), exp(response_to_use)', [holdout_end 1]);

% untransformed response
trunc_ET_mse_test  = fKrigingZ_holdout_MSE_NaNready( trunc_ET_response_est, response_to_use', [holdout_end 1]);
trunc_E_mse_test   = fKrigingZ_holdout_MSE_NaNready( trunc_E_response_est, response_to_use', [holdout_end 1]);
trunc_T_mse_test   = fKrigingZ_holdout_MSE_NaNready( trunc_T_response_est, response_to_use', [holdout_end 1]);
Tikh_ET_mse_test  = fKrigingZ_holdout_MSE_NaNready( Tikh_ET_response_est, response_to_use', [holdout_end 1]);
Tikh_E_mse_test   = fKrigingZ_holdout_MSE_NaNready( Tikh_E_response_est, response_to_use', [holdout_end 1]);
Tikh_T_mse_test   = fKrigingZ_holdout_MSE_NaNready( Tikh_T_response_est, response_to_use', [holdout_end 1]);
ss_tot = fKrigingZ_holdout_MSE_NaNready( response_mean*ones(size(Tikh_ET_response_est)), response_to_use', [holdout_end 1]);


ET_r2_trunc = 1 - trunc_ET_mse_test/ss_tot;
E_r2_trunc = 1 - trunc_E_mse_test/ss_tot;
T_r2_trunc = 1 - trunc_T_mse_test/ss_tot;

ET_r2_Tikh = 1 - Tikh_ET_mse_test/ss_tot;
E_r2_Tikh = 1 - Tikh_E_mse_test/ss_tot;
T_r2_Tikh = 1 - Tikh_T_mse_test/ss_tot;



%
model_name = {"E";"T";"ET"};
trunc_R2 = [E_r2_trunc; T_r2_trunc; ET_r2_trunc];
Tikh_R2 = [E_r2_Tikh; T_r2_Tikh; ET_r2_Tikh];
trunc_MSE = [trunc_E_mse_test; trunc_T_mse_test; trunc_ET_mse_test ];
Tikh_MSE = [Tikh_E_mse_test; Tikh_T_mse_test; Tikh_ET_mse_test ];
table( model_name, trunc_MSE, trunc_R2, Tikh_MSE, Tikh_R2)
