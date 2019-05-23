fileName = 'z_sim_threshold_value.csv';

if (~exist('data','var') )
    opts = detectImportOptions(fileName);
    opts.SelectedVariableNames = {'iCase','seed','transfer_Threshold_value',...
        'transfer_RMSE','ZErrorTrain','ZErrorTrain_trueDynamics_',...
        'ZErrorTest','ZErrorTest_trueDynamics_','ZErrorTrain_holdout_'} ;        
    data = readtable(fileName, opts);
end

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


% my rho
square.nGridTime = [...
    300 300 300 300 300;...
    600 600 600 600 600;...
    900 900 900 900 900];
vec_nGridTime = [300 ;600; 900];
square.nRandi_max = [...
    5 10 20 30 40;...
    5 10 20 30 40;...
    5 10 20 30 40];
my_threshold = 2 .* square.nGridTime.^(-1/3) .* (square.nRandi_max./2).^(-1/2);
my_threshold_disp = [my_threshold  vec_nGridTime.^(-1/2) ];

% read files
for a = 1:3
    for b = 1:5
        
        iCase = (b-1)*3 + a;
        subdata = data( data.iCase == iCase, : );
        
%         find closest rho
        all_thresholds = unique(subdata.transfer_Threshold_value);
        [~,I] = min( abs(my_threshold(a,b) - all_thresholds) );
        my_threshold_this = all_thresholds(I);
        
        subdata = subdata( abs(subdata.transfer_Threshold_value - my_threshold_this) < 0.00001, : );
        
        
        % outputs for the 1st paper        
%         
%         a = find(simulation_cases.nRandi_max(iCase));
%         b = find(simulation_cases.nGridTime(iCase));
        
        square.nSimulations(a,b) = length(unique(subdata.seed));
        square.iCase(a,b) = iCase;
%         
%         square.covLag0RMSE(a,b) = mean( subdata.covLag0RMSE);
%         
%         square.specDensityRMSE_FRO(a,b) = mean( subdata.specDensityRMSE_FRO);
%         square.specDensityRMSE_FRO_sd(a,b) = sqrt(var( subdata.specDensityRMSE_FRO));
%         
%         square.dynamicKrigRMSE(a,b) = median( subdata.dynamicKrigRMSE(subdata.sigma_2EstimLocQv_used_>0.01) );
%         square.dynamicKrigRMSE_IQR(a,b) = iqr( subdata.dynamicKrigRMSE(subdata.sigma_2EstimLocQv_used_>0.01) );
%         %     square.staticKrigRMSE(a,b) = median( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
%         %     square.staticKrigRMSE_IQR(a,b) = iqr( data.staticKrigRMSE(data.sigma_2EstimLocQv_used_>0.05) );
%         
%         square.percentage_good_sigma(a,b) = mean(subdata.sigma_2EstimLocQv_used_>0.01)*100;
%         
%         square.sigma_2EstimLocLin(a,b) = mean( subdata.sigma_2EstimLocLin.^2);
%         square.sigma_2EstimLocQv_used_(a,b) = mean( subdata.sigma_2EstimLocQv_used_.^2);
        
        %     rmse_static  = data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_dynamic = data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05);
        %     rmse_static  = data.staticKrigRMSE;
        %     rmse_dynamic = data.dynamicKrigRMSE;
        %     square.gain(a,b) = median( (rmse_static ./ rmse_dynamic-1)*100 );
        
        % outputs for the 2nd paper
        
        square.transfer_RMSE(a,b) = median(subdata.transfer_RMSE);
        square.ZErrorTrain(a,b) = median(subdata.ZErrorTrain);
        square.ZErrorTrain_trueDynamics_(a,b) = median(subdata.ZErrorTrain_trueDynamics_);
        square.ZErrorTest(a,b) = median(subdata.ZErrorTest);
        square.ZErrorTest_true(a,b) = median(subdata.ZErrorTest_trueDynamics_);
        square.ZErrorTrain_holdout_(a,b) = median(subdata.ZErrorTrain_holdout_);
        
        
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


%%
yvalues = {'300','600','900'};
xvalues = {'2,5','5','10','15','20'};

subplot(3,3,1)
h = heatmap(xvalues,yvalues,square.transfer_RMSE);
h.Title = 'transfer RMSE';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,2)
h = heatmap(xvalues,yvalues,square.ZErrorTrain);
h.Title = 'ZErrorTrain';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,3)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_trueDynamics_);
h.Title = 'ZErrorTrain_trueDynamics_';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,4)
h = heatmap(xvalues,yvalues,square.ZErrorTest);
h.Title = 'ZErrorTest';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,5)
h = heatmap(xvalues,yvalues,square.ZErrorTest_true);
h.Title = 'ZErrorTest_true';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,6)
h = heatmap(xvalues,yvalues,square.ZErrorTrain_holdout_);
h.Title = 'ZErrorTrain_holdout_';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';

subplot(3,3,7)
yvalues2 = {'300','600','900'};
xvalues2 = {'2,5','5','10','15','20','\infty'};
h = heatmap(xvalues2,yvalues2,my_threshold_disp);
h.Title = 'my_threshold';
h.XLabel = 'avg number of obs. per curve';
h.YLabel = 'T';