fileName_Tikhonov = 'z_sim_tikhonov.csv';
fileName_truncation = 'z_sim_threshold_value.csv';
fileName_complete_Tikhonov = 'z_sim_complete_tikhonov.csv';
fileName_complete_truncation = 'z_sim_complete_threshold_multiplier.csv';

dataFull_truncation = readtable(fileName_truncation);
dataFull_Tikhonov = readtable(fileName_Tikhonov);
dataFull_complete_truncation = readtable(fileName_complete_truncation);
dataFull_complete_Tikhonov = readtable(fileName_complete_Tikhonov);

% select simulated process
simulated_process_select = 4;
data_complete_truncation = dataFull_complete_truncation(dataFull_complete_truncation.simulated_process == simulated_process_select, :);
data_complete_Tikhonov = dataFull_complete_Tikhonov(dataFull_complete_Tikhonov.simulated_process == simulated_process_select, :);
data_truncation = dataFull_truncation(dataFull_truncation.simulated_process == simulated_process_select, :);
data_Tikhonov = dataFull_Tikhonov(dataFull_Tikhonov.simulated_process == simulated_process_select, :);


prepare_simulation_cases;



% results identical to the 1st paper
square.covLag0RMSE = nan(3,5);
square.specDensityRMSE_FRO = nan(3,6);
square.specDensityRMSE_FRO_sd = nan(3,6);
square.specCrossDensityRMSE_L2 = nan(3,6);
square.dynamicKrigRMSE = nan(3,6);
square.dynamicKrigRMSE_IQR= nan(3,6);
square.staticKrigRMSE = nan(3,6);
square.staticKrigRMSE_IQR = nan(3,6);
square.sigma_2EstimLocLin = nan(3,6);
square.sigma_2EstimLocQv_used_ = nan(3,6);
square.gain = nan(3,6);
square.percentage_good_sigma = nan(3,6);
square.bw_F = nan(3,6);
square.nRandi_max = nan(3,6);
square.nGridTime = nan(3,6);
square.nSimulations = nan(3,6);
square.iCase = nan(3,6);
% 
% % results for the second paper
% square.transfer_RMSE = nan(3,6);
% square.ZErrorTrain = nan(3,6);
% square.ZErrorTrain_true = nan(3,6);
% square.ZErrorTest = nan(3,6);
% square.ZErrorTest_true = nan(3,6);
% square.ZErrorTrain_holdout_ = nan(3,6);
% 
% square.transfer_RMSE_best_k = nan(3,6);
% square.ZErrorTrain_best_k = nan(3,6);
% square.ZErrorTest_best_k = nan(3,6);
% square.ZErrorTrain_holdout_best_k = nan(3,6);
% square.transfer_RMSE_best_rho = nan(3,6);
% square.ZErrorTrain_best_rho = nan(3,6);
% square.ZErrorTest_best_rho = nan(3,6);
% square.ZErrorTrain_holdout_best_rho = nan(3,6);

% for comparison of truncation and Tikhonov
square_truncation.transfer_RMSE = nan(3,6);
square_truncation.ZErrorTrain = nan(3,6);
square_truncation.ZErrorTrain_true = nan(3,6);
square_truncation.ZErrorTest = nan(3,6);
square_truncation.ZErrorTrain_holdout_ = nan(3,6);
square_truncation.ZErrorTrain_holdoutTrueDynamics_ = nan(3,6);

square_Tikhonov.transfer_RMSE = nan(3,6);
square_Tikhonov.ZErrorTrain = nan(3,6);
square_Tikhonov.ZErrorTrain_true = nan(3,6);
square_Tikhonov.ZErrorTest = nan(3,6);
square_Tikhonov.ZErrorTrain_holdout_ = nan(3,6);
square_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ = nan(3,6);


% read files - sparse
for a = 1:3
    for b = 1:5
        
        iCase = (b-1)*3 + a;
        subdata_truncation = data_truncation( data_truncation.iCase == iCase, : );
        subdata_Tikhonov = data_Tikhonov( data_Tikhonov.iCase == iCase, : );
        
        
        % outputs for the 1st paper        
%         
%         a = find(simulation_cases.nRandi_max(iCase));
%         b = find(simulation_cases.nGridTime(iCase));
        
%         square.nRandi_max(a,b) = subdata_truncation.nRandi_max(1);
%         square.nGridTime(a,b) = subdata_truncation.T(1);
        square.nSimulations(a,b) = length(unique(subdata_truncation.seed));
        square.iCase(a,b) = iCase;
        
        square.covLag0RMSE(a,b) = median( subdata_truncation.covLag0RMSE);
        
        square.specDensityRMSE_FRO(a,b) = median( subdata_truncation.specDensityRMSE_FRO);
        square.specDensityRMSE_FRO_sd(a,b) = sqrt(var( subdata_truncation.specDensityRMSE_FRO));
        
        square.specDensityRMSE_FRO(a,b) = median( subdata_truncation.specDensityRMSE_FRO);
        square.specCrossDensityRMSE_L2(a,b) = meadian( subdata_truncation.specCrossDensityRMSE_L2 );
        
        square.dynamicKrigRMSE(a,b) = median( subdata_truncation.dynamicKrigRMSE(subdata_truncation.sigma_2EstimLocQv_used_>0.01) );
        square.dynamicKrigRMSE_IQR(a,b) = iqr( subdata_truncation.dynamicKrigRMSE(subdata_truncation.sigma_2EstimLocQv_used_>0.01) );
        %     square.staticKrigRMSE(a,b) = median( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        %     square.staticKrigRMSE_IQR(a,b) = iqr( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        
        square.percentage_good_sigma(a,b) = median(subdata_truncation.sigma_2EstimLocQv_used_>0.01)*100;
        
        square.sigma_2EstimLocLin(a,b) = median( subdata_truncation.sigma_2EstimLocLin.^2);
        square.sigma_2EstimLocQv_used_(a,b) = median( subdata_truncation.sigma_2EstimLocQv_used_.^2);
   
        % outputs for the 2nd paper
        square_truncation.transfer_RMSE(a,b) = median( subdata_truncation.transfer_RMSE );
        square_truncation.ZErrorTrain(a,b) = median( subdata_truncation.ZErrorTrain );
        square_truncation.ZErrorTrain_true(a,b) = median( subdata_truncation.ZErrorTrain_trueDynamics_ );
        square_truncation.ZErrorTest(a,b) = median( subdata_truncation.ZErrorTest );
        square_truncation.ZErrorTrain_holdout_(a,b) = median( subdata_truncation.ZErrorTrain_holdout_ );
        square_truncation.ZErrorTrain_holdoutTrueDynamics_(a,b) = median( subdata_truncation.ZErrorTrain_holdoutTrueDynamics_ );
        
        square_Tikhonov.transfer_RMSE(a,b) = median( subdata_Tikhonov.transfer_RMSE );
        square_Tikhonov.ZErrorTrain(a,b) = median( subdata_Tikhonov.ZErrorTrain );
        square_Tikhonov.ZErrorTrain_true(a,b) = median( subdata_Tikhonov.ZErrorTrain_trueDynamics_ );
        square_Tikhonov.ZErrorTest(a,b) = median( subdata_Tikhonov.ZErrorTest );
        square_Tikhonov.ZErrorTrain_holdout_(a,b) = median( subdata_Tikhonov.ZErrorTrain_holdout_ );
        square_Tikhonov.ZErrorTrain_holdoutTrueDynamics_(a,b) = median( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
        
        
    end
end

%% read files - complete
Ts = [300 600 900];
for a = 1:3
    b = 6; % the last column
    nGridTime = Ts(a);
    
    subdata_truncation = data_complete_truncation( data_complete_truncation.T == nGridTime, : );
    subdata_Tikhonov = data_complete_Tikhonov( data_complete_Tikhonov.T == nGridTime, : );
    
    square.nSimulations(a,b) = length(unique(subdata_truncation.seed));
    
    % outputs for the 2nd paper
    square_truncation.transfer_RMSE(a,b) = median( subdata_truncation.transfer_RMSE );
    square_truncation.ZErrorTrain(a,b) = median( subdata_truncation.ZErrorTrain );
    square_truncation.ZErrorTrain_true(a,b) = median( subdata_truncation.ZErrorTrain_trueDynamics_ );
    square_truncation.ZErrorTest(a,b) = median( subdata_truncation.ZErrorTest );
    square_truncation.ZErrorTrain_holdout_(a,b) = median( subdata_truncation.ZErrorTrain_holdout_ );
    square_truncation.ZErrorTrain_holdoutTrueDynamics_(a,b) = median( subdata_truncation.ZErrorTrain_holdoutTrueDynamics_ );
    
    square_Tikhonov.transfer_RMSE(a,b) = median( subdata_Tikhonov.transfer_RMSE );
    square_Tikhonov.ZErrorTrain(a,b) = median( subdata_Tikhonov.ZErrorTrain );
    square_Tikhonov.ZErrorTrain_true(a,b) = median( subdata_Tikhonov.ZErrorTrain_trueDynamics_ );
    square_Tikhonov.ZErrorTest(a,b) = median( subdata_Tikhonov.ZErrorTest );
    square_Tikhonov.ZErrorTrain_holdout_(a,b) = median( subdata_Tikhonov.ZErrorTrain_holdout_ );
    square_Tikhonov.ZErrorTrain_holdoutTrueDynamics_(a,b) = median( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
        
end

% hist( (rmse_static ./ rmse_dynamic-1)*100 )

%% %% create the heat map
square.nRandi_max
square.nGridTime
square.nSimulations
square.iCase

% 
% %%
% yvalues = {'300','600','900'};
% xvalues = {'2,5','5','10','15','20','inf'};
% 
% subplot(3,4,1)
% h = heatmap(xvalues,yvalues,square_truncation.transfer_RMSE);
% h.Title = 'TRUNC: transfer RMSE';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,5)
% h = heatmap(xvalues,yvalues,square_Tikhonov.transfer_RMSE);
% h.Title = 'TIKH: transfer RMSE';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,9)
% h = heatmap(xvalues,yvalues, square_truncation.transfer_RMSE ./ square_Tikhonov.transfer_RMSE );
% h.Title = 'ration TRUNC/TIKH';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% %
% subplot(3,4,2)
% h = heatmap(xvalues,yvalues,square_truncation.ZErrorTest);
% h.Title = 'TRUNC: ZErrorTest';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,6)
% h = heatmap(xvalues,yvalues,square_Tikhonov.ZErrorTest);
% h.Title = 'TIKH: ZErrorTest';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,10)
% h = heatmap(xvalues,yvalues, square_truncation.ZErrorTest ./ square_Tikhonov.ZErrorTest );
% h.Title = 'ration TRUNC/TIKH';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% %
% subplot(3,4,3)
% h = heatmap(xvalues,yvalues,square_truncation.ZErrorTrain_holdout_);
% h.Title = 'TRUNC: ZErrorTrain_holdout_';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,7)
% h = heatmap(xvalues,yvalues,square_Tikhonov.ZErrorTrain_holdout_);
% h.Title = 'TIKH: ZErrorTrain_holdout_';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,11)
% h = heatmap(xvalues,yvalues, square_truncation.ZErrorTrain_holdout_ ./ square_Tikhonov.ZErrorTrain_holdout_ );
% h.Title = 'ration TRUNC/TIKH';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% %
% subplot(3,4,4)
% h = heatmap(xvalues,yvalues,square_truncation.ZErrorTrain_true);
% h.Title = 'TRUNC: ZErrorTrain_true';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,8)
% h = heatmap(xvalues,yvalues,square_Tikhonov.ZErrorTrain_true);
% h.Title = 'TIKH: ZErrorTrain_true';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';
% 
% subplot(3,4,12)
% h = heatmap(xvalues,yvalues, square_truncation.ZErrorTrain_true ./ square_Tikhonov.ZErrorTrain_true );
% h.Title = 'ration TRUNC/TIKH';
% h.XLabel = 'avg number of obs. per curve';
% h.YLabel = 'T';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots fo the paper

data_truncation = dataFull_truncation;
data_Tikhonov = dataFull_Tikhonov;

%% spaghetti plot, dependence on T and N^avg
dir_name = 'figures_paper';
mkdir(dir_name);
h_fig = figure('Position',  [200,200, 1000, 400]);

Ts = [300 600 900];
for ii = 1:3
    T = Ts(ii)';
    subplot(1,3,ii)    
    plot( square_truncation.transfer_RMSE(ii,:) , 'k-x')
    hold on
    plot( square_Tikhonov.transfer_RMSE(ii,:) , 'k--o' )    
    hold off
    ylim([0, 1])
    xlim([0.5 6.5])
    xticks([1 2 3 4 5 6])
    xticklabels({2.5 5 10 15 20 '\infty'})
    xlabel('avg. number of obs. per curve')
    % ylabel('$\delta^{\mathcal{B}}$ - filter coefficients prediction RMSE','Interpreter','latex')
    ylabel('\delta^B - filter coefficients prediction RMSE')
    title(['T = ',num2str(T)])
end

% saveas(h_fig,[dir_name,'/filter1.eps'])
xxx
%% spaghetti plot - prediction, dependence on T and N^avg
Ts = [300 600 900];
for ii = 1:3
    T = Ts(ii)';
    subplot(1,3,ii)    
    plot( square_truncation.ZErrorTrain_holdout_(ii,:) , 'k-x')
    hold on
    plot( square_Tikhonov.ZErrorTrain_holdout_(ii,:) , 'k--o' )
    plot( square_Tikhonov.ZErrorTrain_holdoutTrueDynamics_(ii,:) , 'k:*' )
    hold off
    ylim([0, 1])
    xlim([0.5 6.5])
    xticks([1 2 3 4 5 6])
    xticklabels({2.5 5 10 15 20 '\infty'})
    xlabel('avg. number of obs. per curve')
    % ylabel('$\delta^{\mathcal{B}}$ - filter coefficients prediction RMSE','Interpreter','latex')
    ylabel('\delta^{pred} - filter coefficients prediction RMSE')
    title(['T = ',num2str(T)])
end
% saveas(h_fig,[dir_name,'/pred1.eps'])

%% dependence on the simulated_process
simulated_processes = [0 1 4];
response_processes = [1 2 3];
deltaB_by_sp_truncation = nan(1,3);
deltaB_by_sp_Tikhonov = nan(1,3);
deltaB_by_rp_truncation = nan(1,3);
deltaB_by_rp_Tikhonov = nan(1,3);
deltaPred_by_sp_truncation = nan(1,3);
deltaPred_by_sp_Tikhonov = nan(1,3);
deltaPred_by_sp_true = nan(1,3);
deltaPred_by_rp_truncation = nan(1,3);
deltaPred_by_rp_Tikhonov = nan(1,3);
deltaPred_by_rp_true = nan(1,3);
for ii = 1:3
    % simulated_process
    simulated_process = simulated_processes(ii);
    
    subdata_truncation = data_truncation( data_truncation.simulated_process == simulated_process, : );
	subdata_Tikhonov = data_Tikhonov( data_Tikhonov.simulated_process == simulated_process, : );
    
    deltaB_by_sp_truncation(ii) = median( subdata_truncation.transfer_RMSE );
    deltaB_by_sp_Tikhonov(ii) = median( subdata_Tikhonov.transfer_RMSE );
    deltaPred_by_sp_truncation(ii) = median( subdata_truncation.ZErrorTrain_holdout_ );
    deltaPred_by_sp_Tikhonov(ii) = median( subdata_Tikhonov.ZErrorTrain_holdout_ );
    deltaPred_by_sp_true(ii) = median( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
    
    % response_processes
    response_process = response_processes(ii);
    
    subdata_truncation = data_truncation( data_truncation.response_process == response_process, : );
	subdata_Tikhonov = data_Tikhonov( data_Tikhonov.response_process == response_process, : );
    
    deltaB_by_rp_truncation(ii) = median( subdata_truncation.transfer_RMSE );
    deltaB_by_rp_Tikhonov(ii) = median( subdata_Tikhonov.transfer_RMSE );
    deltaPred_by_rp_truncation(ii) = median( subdata_truncation.ZErrorTrain_holdout_ );
    deltaPred_by_rp_Tikhonov(ii) = median( subdata_Tikhonov.ZErrorTrain_holdout_ );
    deltaPred_by_rp_true(ii) = median( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
end



subplot(1,2,1)
plot( deltaB_by_sp_truncation , 'k-x')
hold on
plot( deltaB_by_sp_Tikhonov , 'k--o' )
hold off
ylim([0, 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'iid','FAR(1)','FMA(4)'})
ylabel('\delta^B - filter coefficients prediction RMSE')
title('Dependence on dynamics of \{X_t\}')

subplot(1,2,2)
plot( deltaB_by_rp_truncation , 'k-x')
hold on
plot( deltaB_by_rp_Tikhonov , 'k--o' )
hold off
ylim([0, 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'(reg1)','(reg2)','(reg3)'})
ylabel('\delta^B - filter coefficients prediction RMSE')
title('Dependence on regression model')

% saveas(h_fig,[dir_name,'/filter2.eps'])

%%
subplot(1,2,1)
plot( deltaPred_by_sp_truncation , 'k-x')
hold on
plot( deltaPred_by_sp_Tikhonov , 'k--o' )
plot( deltaPred_by_sp_true , 'k:*' )
hold off
ylim([0, 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'iid','FAR(1)','FMA(4)'})
ylabel('\delta^{pred} - prediction RMSE')
title('Dependence on dynamics of \{X_t\}')

subplot(1,2,2)
plot( deltaPred_by_rp_truncation , 'k-x')
hold on
plot( deltaPred_by_rp_Tikhonov , 'k--o' )
plot( deltaPred_by_rp_true , 'k:*' )
hold off
ylim([0, 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'(reg1)','(reg2)','(reg3)'})
ylabel('\delta^{pred} - prediction RMSE')
title('Dependence on regression model')

% saveas(h_fig,[dir_name,'/pred2.eps'])




