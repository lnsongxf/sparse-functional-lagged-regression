fileName = 'z_sim_threshold_multiplier.csv';

 data = readtable(fileName);

prepare_simulation_cases;

% results identical to the 1st paper
square.covLag0RMSE = nan(3,5);
square.specDensityRMSE_FRO = nan(3,5);
square.specDensityRMSE_FRO_sd = nan(3,5);
square.dynamicKrigRMSE = nan(3,5);
square.dynamicKrigRMSE_IQR= nan(3,5);
square.staticKrigRMSE = nan(3,5);
square.staticKrigRMSE_IQR = nan(3,5);
square.sigma_2EstimLocLin = nan(3,5);
square.sigma_2EstimLocQv_used_ = nan(3,5);
square.gain = nan(3,5);
square.percentage_good_sigma = nan(3,5);
square.bw_F = nan(3,5);
square.nRandi_max = nan(3,5);
square.nGridTime = nan(3,5);
square.nSimulations = nan(3,5);
square.iCase = nan(3,5);

% results for the second paper
square.transfer_RMSE = nan(3,5);
square.ZErrorTrain = nan(3,5);
square.ZErrorTrain_true = nan(3,5);
square.ZErrorTest = nan(3,5);
square.ZErrorTest_true = nan(3,5);
square.ZErrorTrain_holdout_ = nan(3,5);

square.transfer_RMSE_best_k = nan(3,5);
square.ZErrorTrain_best_k = nan(3,5);
square.ZErrorTest_best_k = nan(3,5);
square.ZErrorTrain_holdout_best_k = nan(3,5);
square.transfer_RMSE_best_rho = nan(3,5);
square.ZErrorTrain_best_rho = nan(3,5);
square.ZErrorTest_best_rho = nan(3,5);
square.ZErrorTrain_holdout_best_rho = nan(3,5);

% read files
for a = 1:3
    for b = 1:5
        
        iCase = (b-1)*3 + a;
        subdata = data( data.iCase == iCase, : );
        
        
        % outputs for the 1st paper        
%         
%         a = find(simulation_cases.nRandi_max(iCase));
%         b = find(simulation_cases.nGridTime(iCase));
        
        square.nRandi_max(a,b) = subdata.nRandi_max(1);
        square.nGridTime(a,b) = subdata.T(1);
        square.nSimulations(a,b) = length(unique(subdata.seed));
        square.iCase(a,b) = iCase;
        
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
        multipliers = unique(subdata.transfer_Threshold_multiplier);
        multipliers_n = length(multipliers);
                
        vec = [];
        vec.transfer_RMSE = zeros(1,multipliers_n);
        vec.ZErrorTrain = zeros(1,multipliers_n);
        vec.ZErrorTrain_true = zeros(1,multipliers_n);
        vec.ZErrorTest = zeros(1,multipliers_n);
        vec.ZErrorTest_true = zeros(1,multipliers_n);
        vec.ZErrorTrain_holdout_ = zeros(1,multipliers_n);
        
        for k = 1:multipliers_n
            rho = multipliers(k);
            vec.transfer_RMSE(k) = median( subdata.transfer_RMSE( subdata.transfer_Threshold_multiplier == rho ) );
            vec.ZErrorTrain(k) = median( subdata.ZErrorTrain( subdata.transfer_Threshold_multiplier == rho ) );
            vec.ZErrorTrain_true(k) = median( subdata.ZErrorTrain_trueDynamics_( subdata.transfer_Threshold_multiplier == rho ) );
            vec.ZErrorTest(k) = median( subdata.ZErrorTest( subdata.transfer_Threshold_multiplier == rho ) );
            vec.ZErrorTest_true(k) = median( subdata.ZErrorTest_trueDynamics_( subdata.transfer_Threshold_multiplier == rho ) );
            vec.vec.ZErrorTrain_holdout_(k) = median( subdata.ZErrorTrain_holdout_( subdata.transfer_Threshold_multiplier == rho ) );
        end
        
        [square.transfer_RMSE(a,b), square.transfer_RMSE_best_k(a,b)] = min(vec.transfer_RMSE);
        [square.ZErrorTrain(a,b), square.ZErrorTrain_best_k(a,b)] = min(vec.ZErrorTrain);
        [square.ZErrorTrain_true(a,b), ~] = min(vec.ZErrorTrain_true);
        [square.ZErrorTest(a,b), square.ZErrorTest_best_k(a,b)] = min(vec.ZErrorTest);
        [square.ZErrorTest_true(a,b), ~] = min(vec.ZErrorTest_true);
        [square.ZErrorTrain_holdout_(a,b), square.ZErrorTrain_holdout_best_k(a,b)] = min(vec.ZErrorTest);
        
        square.transfer_RMSE_best_multiplier(a,b) = multipliers(square.transfer_RMSE_best_k(a,b));
        square.ZErrorTrain_best_multiplier(a,b) = multipliers(square.ZErrorTrain_best_k(a,b));
        square.ZErrorTest_best_multiplier(a,b) = multipliers(square.ZErrorTest_best_k(a,b));
        square.ZErrorTrain_holdout_best_multiplier(a,b) = multipliers(square.ZErrorTrain_holdout_best_k(a,b));
        
        
    end
end

% hist( (rmse_static ./ rmse_dynamic-1)*100 )

%% %% create the heat map
square.nRandi_max
square.nGridTime
square.nSimulations
square.iCase


%%
yvalues = {'300','600','900'};
xvalues = {'2,5','5','10','15','20'};

subplot(3,3,1)
h = heatmap(xvalues,yvalues,square.transfer_RMSE);
h.Title = 'transfer RMSE';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,2)
h = heatmap(xvalues,yvalues,square.transfer_RMSE_best_multiplier);
h.Title = 'transfer RMSE - best multiplier';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,4)
h = heatmap(xvalues,yvalues,square.ZErrorTrain);
h.Title = 'ZErrorTrain';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,7)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_best_multiplier);
h.Title = 'ZErrorTrain_best_multiplier';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,5)
h = heatmap(xvalues,yvalues,square.ZErrorTest);
h.Title = 'ZErrorTest';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,8)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_best_multiplier);
h.Title = 'ZErrorTest_best_multiplier';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,6)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_holdout_);
h.Title = 'ZErrorTrain_holdout_';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';


subplot(3,3,9)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_holdout_best_multiplier);
h.Title = 'ZErrorTrain_holdout_best_multiplier';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,3)
h = heatmap(xvalues,yvalues,square.nSimulations);
h.Title = 'nSimulations';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';



