fileName = 'z_sim_all.csv';

% data = readtable(fileName);

prepare_simulation_cases;

% results identical to the 1st paper
square.covLag0RMSE = nan(6,9);
square.specDensityRMSE_FRO = nan(6,9);
square.specDensityRMSE_FRO_sd = nan(6,9);
square.dynamicKrigRMSE = nan(6,9);
square.dynamicKrigRMSE_IQR= nan(6,9);
square.staticKrigRMSE = nan(6,9);
square.staticKrigRMSE_IQR = nan(6,9);
square.sigma_2EstimLocLin = nan(6,9);
square.sigma_2EstimLocQv_used_ = nan(6,9);
square.gain = nan(6,9);
square.percentage_good_sigma = nan(6,9);
square.bw_F = nan(6,9);
square.nRandi_max = nan(6,9);
square.nGridTime = nan(6,9);
square.nSimulations = nan(6,9);
square.iCase = nan(6,9);

% results for the second paper
square.transfer_RMSE = nan(6,9);
square.ZErrorTrain = nan(6,9);
square.ZErrorTrain_true = nan(6,9);
square.ZErrorTest = nan(6,9);
square.ZErrorTest_true = nan(6,9);

square.transfer_RMSE_best_k = nan(6,9);
square.ZErrorTrain_best_k = nan(6,9);
square.ZErrorTest_best_k = nan(6,9);

simulated_processes = 0:5;
response_processes = 1:9;

% read files
for simulated_process = simulated_processes
    for response_process = response_processes
        
        a = simulated_process+1;
        b = response_process;
        
        subdata = data( data.simulated_process == simulated_process & data.response_process==response_process, : );
                
        square.nSimulations(a,b) = length(unique(subdata.seed));        
        square.covLag0RMSE(a,b) = mean( subdata.covLag0RMSE);        
        square.specDensityRMSE_FRO(a,b) = mean( subdata.specDensityRMSE_FRO);
        square.specDensityRMSE_FRO_sd(a,b) = sqrt(var( subdata.specDensityRMSE_FRO));
        
        square.dynamicKrigRMSE(a,b) = median( subdata.dynamicKrigRMSE(subdata.sigma_2EstimLocQv_used_>0.01) );
        square.dynamicKrigRMSE_IQR(a,b) = iqr( subdata.dynamicKrigRMSE(subdata.sigma_2EstimLocQv_used_>0.01) );
        %     square.staticKrigRMSE(a,b) = median( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        %     square.staticKrigRMSE_IQR(a,b) = iqr( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        
        square.percentage_good_sigma(a,b) = mean(subdata.sigma_2EstimLocQv_used_>0.01)*100;
        
        square.sigma_2EstimLocLin(a,b) = mean( subdata.sigma_2EstimLocLin.^2);
        square.sigma_2EstimLocQv_used_(a,b) = mean( subdata.sigma_2EstimLocQv_used_.^2);
        
        %     rmse_static  = data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_dynamic = data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_static  = data.staticKrigRMSE;
        %     rmse_dynamic = data.dynamicKrigRMSE;
        %     square.gain(a,b) = median( (rmse_static ./ rmse_dynamic-1)*100 );
        
        % outputs for the 2nd paper
        trunc_max = max(subdata.transfer_TruncationLevel_nEigenvalues_);
        
        vec = [];
        vec.transfer_RMSE = zeros(1,trunc_max);
        vec.ZErrorTrain = zeros(1,trunc_max);
        vec.ZErrorTrain_true = zeros(1,trunc_max);
        vec.ZErrorTest = zeros(1,trunc_max);
        vec.ZErrorTest_true = zeros(1,trunc_max);
        
        for k = 1:trunc_max
            vec.transfer_RMSE(k) = median( subdata.transfer_RMSE( subdata.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTrain(k) = median( subdata.ZErrorTrain( subdata.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTrain_true(k) = median( subdata.ZErrorTrain_trueDynamics_( subdata.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTest(k) = median( subdata.ZErrorTest( subdata.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTest_true(k) = median( subdata.ZErrorTest_trueDynamics_( subdata.transfer_TruncationLevel_nEigenvalues_ == k ) );
        end
        
        [square.transfer_RMSE(a,b), square.transfer_RMSE_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTrain(a,b), square.ZErrorTrain_best_k(a,b)] = min(vec.ZErrorTrain);
        [square.ZErrorTrain_true(a,b), ~] = min(vec.ZErrorTrain_true);
        [square.ZErrorTest(a,b), square.ZErrorTest_best_k(a,b)] = min(vec.ZErrorTest);
        [square.ZErrorTest_true(a,b), ~] = min(vec.ZErrorTest_true);
        
    end
end

% hist( (rmse_static ./ rmse_dynamic-1)*100 )

%% %% create the heat map
square.nSimulations


%%
yvalues = simulated_processes;
xvalues = response_processes;

subplot(3,2,1)
h = heatmap(xvalues,yvalues,square.transfer_RMSE);
h.Title = 'transfer RMSE';
h.XLabel = 'response process';
h.YLabel = 'simulated process';


subplot(3,2,3)
h = heatmap(xvalues,yvalues,square.ZErrorTrain);
h.Title = 'ZErrorTrain';
h.XLabel = 'response process';
h.YLabel = 'simulated process';

subplot(3,2,4)
h = heatmap(xvalues,yvalues,square.ZErrorTest);
h.Title = 'ZErrorTest';
h.XLabel = 'response process';
h.YLabel = 'simulated process';

subplot(3,2,5)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_true);
h.Title = 'ZErrorTrain_true';
h.XLabel = 'response process';
h.YLabel = 'simulated process';

subplot(3,2,6)
h = heatmap(xvalues,yvalues,square.ZErrorTest_true);
h.Title = 'ZErrorTest_true';
h.XLabel = 'response process';
h.YLabel = 'simulated process';

subplot(3,2,2)
h = heatmap(xvalues,yvalues,square.nSimulations);
h.Title = 'nSimulations';
h.XLabel = 'response process';
h.YLabel = 'simulated process';

% %%
% yvalues = simulated_processes;
% xvalues = response_processes;
% 
% subplot(3,4,1)
% h = heatmap(xvalues,yvalues,square.transfer_RMSE);
% h.Title = 'transfer RMSE';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,2)
% h = heatmap(xvalues,yvalues,square.transfer_RMSE_best_k);
% h.Title = 'transfer RMSE - best k';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,5)
% h = heatmap(xvalues,yvalues,square.ZErrorTrain);
% h.Title = 'ZErrorTrain';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,6)
% h = heatmap(xvalues,yvalues,square.ZErrorTrain_best_k);
% h.Title = 'ZErrorTrain_best_k';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,7)
% h = heatmap(xvalues,yvalues,square.ZErrorTest);
% h.Title = 'ZErrorTest';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,8)
% h = heatmap(xvalues,yvalues,square.ZErrorTrain_best_k);
% h.Title = 'ZErrorTest_best_k';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,9)
% h = heatmap(xvalues,yvalues,square.ZErrorTrain_true);
% h.Title = 'ZErrorTrain_true';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% 
% subplot(3,4,11)
% h = heatmap(xvalues,yvalues,square.ZErrorTest_true);
% h.Title = 'ZErrorTest_true';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';
% 
% subplot(3,4,3)
% h = heatmap(xvalues,yvalues,square.nSimulations);
% h.Title = 'nSimulations';
% h.XLabel = 'response process';
% h.YLabel = 'simulated process';



