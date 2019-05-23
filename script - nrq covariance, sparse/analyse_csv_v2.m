fileNamePrefix = 'sim_91_case_';
fileNameSuffix = '.csv';

% first and second matrix to display
square_covLag0RMSE = nan(5,6);
square_specDensityRMSE_FRO = nan(5,6);
square_specDensityRMSE_FRO_sd = nan(5,6);
square_dynamicKrigRMSE = nan(5,6);
square_dynamicKrigRMSE_IQR= nan(5,6);
square_staticKrigRMSE = nan(5,6);
square_staticKrigRMSE_IQR = nan(5,6);
square_sigma_estLocLin = nan(5,6);
square_sigma_estLocQv = nan(5,6);
square_sigma2_estLocLin = nan(5,6);
square_sigma2_estLocQv = nan(5,6);
square_gain = nan(5,6);
square_percentage_good_sigma = nan(5,6);

square_bw_F = nan(5,6);

square_nRandi_max = nan(5,6);
square_nGridTime = nan(5,6);
square_nSimulations = nan(5,6);

% map iCase to the corresponding nGridTime and nRandi_max
simulation_variables.nGridTime = [150 300 450 600 900 1200];
simulation_variables.nRandi_max = [5 10 20 30 40];
nCases = length(simulation_variables.nGridTime) * length(simulation_variables.nRandi_max);
simulation_cases = [];
simulation_cases.nGridTime = zeros(1,nCases);
simulation_cases.nRandi_max = zeros(1,nCases);

iCase = 1;
for i1 = 1:length(simulation_variables.nRandi_max)
    for i2 = 1:length(simulation_variables.nGridTime)
        simulation_cases.nRandi_max_i(iCase) = i1;    
        simulation_cases.nRandi_max_val(iCase) = simulation_variables.nRandi_max(i1);
        simulation_cases.nGridTime_i(iCase) = i2;            
        simulation_cases.nGridTime_i_val(iCase) = simulation_variables.nGridTime(i2);            
        iCase = iCase + 1;
    end
end

% save values for regression
spd.fro = [];
spd.intmax = [];
spd.T = [];
spd.N = [];


% read files
for iCase = 1:30
    fileName = [fileNamePrefix, num2str(iCase), fileNameSuffix];
    data = readtable(fileName);
           
    a = simulation_cases.nRandi_max_i(iCase);
    b = simulation_cases.nGridTime_i(iCase);
    
    square_nRandi_max(a,b) = data.nRandi_max(1);
    square_nGridTime(a,b) = data.T(1);
    square_nSimulations(a,b) = length(unique(data.seed));
    square_iCase(a,b) = iCase;
    
    square_covLag0RMSE(a,b) = mean( data.covLag0RMSE);
    
    square_specDensityRMSE_FRO(a,b) = mean( data.specDensityRMSE_FRO);
    square_specDensityRMSE_FRO_sd(a,b) = sqrt(var( data.specDensityRMSE_FRO));
    
    square_dynamicKrigRMSE(a,b) = median( data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05) );
    square_dynamicKrigRMSE_IQR(a,b) = iqr( data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05) );
    square_staticKrigRMSE(a,b) = median( data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05) );    
    square_staticKrigRMSE_IQR(a,b) = iqr( data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05) );  
    
    square_percentage_good_sigma(a,b) = mean(data.sigma_estLocQv_used_>0.05)*100;
    
    square_sigma_estLocLin(a,b) = mean( data.sigma_estLocLin);
    square_sigma_estLocQv(a,b) = mean( data.sigma_estLocQv_used_);
    square_sigma2_estLocLin(a,b) = mean( data.sigma_estLocLin.^2);
    square_sigma2_estLocQv(a,b) = mean( data.sigma_estLocQv_used_.^2);
    
    rmse_static  = data.staticKrigRMSE(data.sigma_estLocQv_used_>0.05);
    rmse_dynamic = data.dynamicKrigRMSE(data.sigma_estLocQv_used_>0.05);    
%     rmse_static  = data.staticKrigRMSE;
%     rmse_dynamic = data.dynamicKrigRMSE;
    
    square_gain(a,b) = median( (rmse_static ./ rmse_dynamic-1)*100 );
    
%         
%     % data for linear regression
%     spd.fro = [spd.fro, data.specDensityRMSE_FRO( data.q_bartlett == q )' ];
%     spd.intmnax = [spd.fro, data.specDensityRMSE_Intmax( data.q_bartlett == q )' ];
%     spd.T = [spd.T, data.T( data.q_bartlett == q )' ];
%     spd.N = [spd.N, data.avgObsPerCurve( data.q_bartlett == q )' ];
end

% hist( (rmse_static ./ rmse_dynamic-1)*100 )

%% %% create the heat map
square_nRandi_max
square_nGridTime
square_nSimulations
square_iCase


%%

X1 = square_nRandi_max(:)/2;
X2 = square_nGridTime(:);
Y = square_specDensityRMSE_FRO(:)


Xlog = [ones(length(X1),1) log(X1) log(X2)];
Xlogf = Xlog;
Ylogf = log(Y);
% Xlogf = Xlog([1:24 26:29],:);
% Ylogf = log(Y([1:24 26:29],:));
b = Xlogf\Ylogf % the last coefficient (with nGridTime) should be smaller

Ylogf_fit = Xlogf*b;
Rsq = 1 - sum((Ylogf - Ylogf_fit).^2)/sum((Ylogf - mean(Ylogf)).^2)


%% log plot
h_fig = figure('Position', [200,200, 800, 600]);


x1_plot = linspace(0.8,3.2,100);
x2_plot = linspace(4.5,7.5,100);
y_plot =  b(1) + b(2)*x1_plot + b(3)*x2_plot' ;
scatter3(   Xlogf(:,2),  Xlogf(:,3)  , Ylogf, 80,  'r', 'filled' )
view(45+90-20, 30)
text(Xlogf(:,2) + 0.08,  Xlogf(:,3)  , Ylogf, num2str( round(exp(Ylogf),2) ),...
    'VerticalAlignment','bottom','HorizontalAlignment','right')
hold on
% scatter3(   Xlogf(:,2),  Xlogf(:,3)  , Ylogf_fit, 100,  'xr' )
for ii = 1:length(Ylogf_fit)
    line( [Xlogf(ii,2),Xlogf(ii,2)], [Xlogf(ii,3),Xlogf(ii,3)], [-1e6,Ylogf(ii)], 'Color', 'k','LineStyle','--', 'LineWidth',0.1 )
    line( [Xlogf(ii,2),Xlogf(ii,2)], [Xlogf(ii,3),Xlogf(ii,3)], [Ylogf_fit(ii),Ylogf(ii)], 'Color', 'k', 'LineWidth',3 )
end
surf( x1_plot, x2_plot, y_plot, 'FaceAlpha',0.5 , 'EdgeColor', 'none', 'FaceColor', 'b')
hold off

title('Spectral density estimation of AR(1), c=0.7') % XXX

ylabel('Time-series length T (log scale)')
yticks( log([150 300 450 600 900 1200]) )
yticklabels([150 300 450 600 900 1200])
ylim(log([120 1300]))

xlabel('N^{max} (log scale)')
xticks( log([2.5 5 10 15 20]) )
xticklabels([5 10 20 30 40])
xlim(log([2 22]))

zlabel('Spec. density estimator RMSE (log scale)')
z_vector = [0.001 0.005 0.01 0.05 0.06 0.07 0.08 0.09 0.1  0.15 0.2 0.3 0.4 0.5 ];
zticks( log(z_vector) )
zticklabels(z_vector)
zlim([-3 0])

set(h_fig,'Units','Inches');
pos = get(h_fig,'Position');
set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_fig, 'AR1_07_spec_density.pdf','-dpdf','-r0') % XXX
close(h_fig)

%%
subplot(3,3,9)

yvalues = {'2.5','5','10','15','20'};
xvalues = {'150','300','450','600','900','1200'};
h = heatmap(xvalues,yvalues,square_covLag0RMSE);

h.Title = 'cov lag 0 RMSE';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,6)
h = heatmap(xvalues,yvalues,square_specDensityRMSE_FRO);
h.Title = 'spec density RMSE (frobenius)';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,1)
h = heatmap(xvalues,yvalues,real(square_dynamicKrigRMSE));
h.Title = 'dynamic kriging RMSE';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,2)
h = heatmap(xvalues,yvalues,real(square_staticKrigRMSE));
h.Title = 'static kriging RMSE';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,3) % gain
% gain = (real(square_staticKrigRMSE) ./ real(square_dynamicKrigRMSE) - 1)*100;
h = heatmap(xvalues,yvalues, square_gain );
h.Title = 'gain';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';


subplot(3,3,7)
h = heatmap(xvalues,yvalues,square_sigma_estLocLin);
h.Title = 'sigma est (loc-lin)';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,8)
h = heatmap(xvalues,yvalues,square_sigma_estLocQv);
h.Title = 'sigma est (loc-qv)';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,4)
h = heatmap(xvalues,yvalues,square_sigma2_estLocLin);
h.Title = 'sigma^2 est (loc-lin)';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';

subplot(3,3,5)
h = heatmap(xvalues,yvalues,square_sigma2_estLocQv);
h.Title = 'sigma^2 est (loc-qv)';
h.YLabel = 'avg number of obs. per curve';
h.XLabel = 'T';



sigma_true = data.sigma(1)
sigma2_true = data.sigma(1)^2

suptitle('FAR(1), c=0.7') % XXX

square_specDensityRMSE_FRO
square_specDensityRMSE_FRO_sd


%% table for kriging


[square_dynamicKrigRMSE(:), square_dynamicKrigRMSE_IQR(:),...
    square_staticKrigRMSE(:), square_staticKrigRMSE_IQR(:),...
    square_gain(:), square_percentage_good_sigma(:)]
