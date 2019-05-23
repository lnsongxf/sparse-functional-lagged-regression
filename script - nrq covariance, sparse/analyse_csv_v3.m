prepare_simulation_cases;

simulated_process = 0; % XXX
response_process = 1; % XXX




% results identical to the 1st paper
square.covLag0RMSE = nan(3,4);
square.specDensityRMSE_FRO = nan(3,4);
square.specDensityRMSE_FRO_sd = nan(3,4);
square.dynamicKrigRMSE = nan(3,4);
square.dynamicKrigRMSE_IQR= nan(3,4);
square.staticKrigRMSE = nan(3,4);
square.staticKrigRMSE_IQR = nan(3,4);
square.sigma_2EstimLocLin = nan(3,4);
square.sigma_2EstimLocQv_used_ = nan(3,4);
square.gain = nan(3,4);
square.percentage_good_sigma = nan(3,4);
square.bw_F = nan(3,4);
square.nRandi_max = nan(3,4);
square.nGridTime = nan(3,4);
square.nSimulations = nan(3,4);
square.iCase = nan(3,4);

% results for the second paper
square.transfer_RMSE = nan(3,4);
square.ZErrorTrain = nan(3,4);
square.ZErrorTrain_true = nan(3,4);
square.ZErrorTest = nan(3,4);
square.ZErrorTest_true = nan(3,4);

square.transfer_RMSE_best_k = nan(3,4);
square.ZErrorTrain_best_k = nan(3,4);
square.ZErrorTrain_true_best_k = nan(3,4);
square.ZErrorTest_best_k = nan(3,4);
square.ZErrorTest_true_best_k = nan(3,4);

% read files
for a = 1:3
    for b = 1:4
        
        iCase = (b-1)*3 + a;
        
        % outputs for the 1st paper
        fileName = ['z_sim_01_c', num2str(iCase,'%02d'),'_p',num2str(simulated_process),'_r',num2str(response_process,'%02d'), '.csv'];
        data = readtable(fileName);
%         
%         a = find(simulation_cases.nRandi_max(iCase));
%         b = find(simulation_cases.nGridTime(iCase));
        
        square.nRandi_max(a,b) = data.nRandi_max(1);
        square.nGridTime(a,b) = data.T(1);
        square.nSimulations(a,b) = length(unique(data.seed));
        square.iCase(a,b) = iCase;
        
        square.covLag0RMSE(a,b) = mean( data.covLag0RMSE);
        
        square.specDensityRMSE_FRO(a,b) = mean( data.specDensityRMSE_FRO);
        square.specDensityRMSE_FRO_sd(a,b) = sqrt(var( data.specDensityRMSE_FRO));
        
        square.dynamicKrigRMSE(a,b) = median( data.dynamicKrigRMSE(data.sigma_2EstimLocQv_used_>0.01) );
        square.dynamicKrigRMSE_IQR(a,b) = iqr( data.dynamicKrigRMSE(data.sigma_2EstimLocQv_used_>0.01) );
        %     square.staticKrigRMSE(a,b) = median( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        %     square.staticKrigRMSE_IQR(a,b) = iqr( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
        
        square.percentage_good_sigma(a,b) = mean(data.sigma_2EstimLocQv_used_>0.01)*100;
        
        square.sigma_2EstimLocLin(a,b) = mean( data.sigma_2EstimLocLin.^2);
        square.sigma_2EstimLocQv_used_(a,b) = mean( data.sigma_2EstimLocQv_used_.^2);
        
        %     rmse_static  = data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_dynamic = data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_static  = data.staticKrigRMSE;
        %     rmse_dynamic = data.dynamicKrigRMSE;
        %     square.gain(a,b) = median( (rmse_static ./ rmse_dynamic-1)*100 );
        
        % outputs for the 2nd paper
        trunc_max = max(data.transfer_TruncationLevel_nEigenvalues_);
        
        vec = [];
        vec.transfer_RMSE = zeros(1,trunc_max);
        vec.ZErrorTrain = zeros(1,trunc_max);
        vec.ZErrorTrain_true = zeros(1,trunc_max);
        vec.ZErrorTest = zeros(1,trunc_max);
        vec.ZErrorTest_true = zeros(1,trunc_max);
        
        for k = 1:trunc_max
            vec.transfer_RMSE(k) = mean( data.transfer_RMSE( data.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTrain(k) = mean( data.ZErrorTrain( data.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTrain_true(k) = mean( data.ZErrorTrain_trueDynamics_( data.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTest(k) = mean( data.ZErrorTest( data.transfer_TruncationLevel_nEigenvalues_ == k ) );
            vec.ZErrorTest_true(k) = mean( data.ZErrorTest_trueDynamics_( data.transfer_TruncationLevel_nEigenvalues_ == k ) );
        end
        
        [square.transfer_RMSE(a,b), square.transfer_RMSE_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTrain(a,b), square.ZErrorTrain_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTrain_true(a,b), square.ZErrorTrain_true_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTest(a,b), square.ZErrorTest_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTest_true(a,b), square.ZErrorTest_true_best_k(a,b)] = min(vec.transfer_RMSE);
        
    end
end

% hist( (rmse_static ./ rmse_dynamic-1)*100 )

%% %% create the heat map
square.nRandi_max
square.nGridTime
square.nSimulations
square.iCase


%%
yvalues = {'250','500','1000'};
xvalues = {'5','10','20','30'};

subplot(3,3,1)
h = heatmap(xvalues,yvalues,square.transfer_RMSE);
h.Title = 'transfer RMSE';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,2)
h = heatmap(xvalues,yvalues,square.transfer_RMSE_best_k);
h.Title = 'transfer RMSE - best k';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,3)
h = heatmap(xvalues,yvalues,square.ZErrorTrain);
h.Title = 'ZErrorTrain';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,4)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_best_k);
h.Title = 'ZErrorTrain_best_k';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,5)
h = heatmap(xvalues,yvalues,square.ZErrorTest);
h.Title = 'ZErrorTest';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,6)
h = heatmap(xvalues,yvalues,square.ZErrorTest_best_k);
h.Title = 'ZErrorTest_best_k';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,7)
h = heatmap(xvalues,yvalues,square.ZErrorTest_true);
h.Title = 'ZErrorTest_true';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,8)
h = heatmap(xvalues,yvalues,square.ZErrorTest_true_best_k);
h.Title = 'ZErrorTest_true_best_k';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,9)
h = heatmap(xvalues,yvalues,square.nSimulations);
h.Title = 'nSimulations';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';
