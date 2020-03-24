%{

This script loads the data "Wank_data.csv" and performs the spectral
analysis of the regressor time series of atmospheric electricity.
- it fits the spectral density from sparse observations
- performs the functional data recovery

WARNING: the runtime of this script requires a lot of memory. We performed
it on a cluster computer.

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}
clear all

addpath('../master')
addpath('../master/fdaM')

show_plots = 0;

dataFull = readtable('Wank_data.csv');


% number of days
nGridTime = length(unique(dataFull.year_day));

% extract the days
Eat_daily = nan(nGridTime, 24);
Eat_indic = nan(nGridTime, 24);
vis_daily = nan(nGridTime, 24);


day = 0;
DOY_last = -9999;
for t = 1:height(dataFull)    
    
    DOY = dataFull.DOY(t)+1;    
    if DOY ~= DOY_last, day = day+1; end    
    DOY_last = DOY;
    hour = dataFull.hour(t)+1;    
    
    Eat_daily(day,hour) = dataFull.E_V_m_Wank(t);    
    vis_daily(day,hour) = dataFull.Vis_km_Wank(t);
    
    % indicator
    Eat_indic(day,hour) = dataFull.E_V_m_Wank_valid_obs(t);    
end

% replace all NaN in the indicator with zero
Eat_indic(isnan(Eat_indic))=0;




%% create ONB
% nGridSpace = 241, nKnots = 24
onb = fPrepareONB(241, 24);



%% prepare the censor dataset
censor = [];
censor.nGridTime = nGridTime;
censor.onb.nGridSpace = onb.nGridSpace;
censor.onb.gridSpace = onb.gridSpace;

% grid for 24 dimensional vector
censor.onb.grid24 = ((1:24)-1)*10 + 6;

censor.zeroone = zeros(nGridTime,onb.nGridSpace); % censorBool(t,x) = 1 iff the I do observe position "x" on the "t"-th curve
censor.Hspace = cell(nGridTime,1);
censor.Honb = cell(nGridTime,1);
censor.grid = cell(nGridTime,1);
censor.nGrid = zeros(nGridTime,1);
censor.data = cell(nGridTime,1);
for t = 1:nGridTime
    censor.nGrid(t,1) = sum( Eat_indic(t,:) );
    for ii = 1:24
        if (Eat_indic(t,ii) == 1)
            % ii = hour, 1 -> 6, 2-> 16 ..... point = (hour-1)*10 + 6
            % inverse conversion: 
            censor.zeroone(t, (ii-1)*10 + 6 ) = 1;
        end
    end
    
    % (actual) current censor matrix H
    numObs = censor.nGrid(t,1);
    Hspace_act = zeros(numObs, onb.nGridSpace );
    grid_act = zeros(numObs,1);
    j = 1;
    for x=1:onb.nGridSpace        
        if (censor.zeroone(t,x) == 1)
            Hspace_act(j,x) = 1;
            grid_act(j) = x;
            j = j+1;
        end        
    end    
    
    censor.Hspace{t,:} = Hspace_act;
    censor.Honb{t,:} = Hspace_act*onb.onbMatrix;
    censor.grid{t,:} = grid_act;
    censor.data{t,:} = Eat_daily(t, (grid_act-6)/10 + 1 );
end

%% info about Eat_daily

sum(censor.nGrid) / nGridTime

%% prepare the scalar response process

response_dailymean = nan( nGridTime, 1 );
response_onehour = nan( nGridTime, 1 );

for t = 1:nGridTime
    
    x = vis_daily(t,:);    
    
    if sum(isnan(x)) < 24
        response_dailymean(t) = mean( x, 'omitnan' );
    end
    
    response_onehour(t) = x(15);
    
end

response_dailymean_available = ~isnan(response_dailymean);
response_onehour_available = ~isnan(response_onehour);

response = response_dailymean;
response_available = response_dailymean_available;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% investigate correlation


Eat_daily_mean = nan(nGridTime, 1);
for t = 1:nGridTime
    Eat_daily_mean(t) = mean( Eat_daily(t,:) ,'omitnan' );
end

if show_plots
    plot(Eat_daily_mean)
    hold on
    plot(response_dailymean)
    hold off
end


%%

if show_plots
    scatter( Eat_daily_mean, response ) 
    m = cov( Eat_daily_mean, response, 'omitrows' );
    c = m(1,2) / sqrt(m(1,1)) / sqrt(m(2,2))
end

%% some plots
% nDays = 5*8;
% offset = 40 * 5;
% for ii = 1:nDays
%     t = offset+ii;
%     subplot(5,8,ii)
%     plot( Eat_daily(t,:) )    
%     hold on
%     nonzero_positions = find(Eat_indic(t,:));
%     scatter( nonzero_positions, Eat_daily(t, nonzero_positions), 'xb' )
%     plot( vis_daily(t,:) )
%     hold off
%     title("day="+num2str(t))
% end


%% some plots - censor

if show_plots
    nDays = 5*8;
    offset = 40 * 11;
    for ii = 1:nDays
        t = offset+ii;
        subplot(5,8,ii)
        plot( censor.onb.gridSpace(censor.onb.grid24), Eat_daily(t,:) )    
        hold on
        nonzero_positions = find(Eat_indic(t,:));
        scatter( censor.onb.gridSpace(censor.grid{t}) , censor.data{t,:}, 'xb' )
        plot( censor.onb.gridSpace(censor.onb.grid24), vis_daily(t,:) )
        if response_available(t)
            plot( censor.onb.gridSpace, response(t) * ones(size(censor.onb.gridSpace)))
        end    
        hold off
        ylim([-50 300])
        title("day="+num2str(t))
    end
end



%% periodicity
scatter( dataFull.DOY, dataFull.E_V_m_Wank )
refline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare the partition for the K-FOLD cross-validation

% how many K batches
censor.cv_K = 10;

% NEW - assign the same CV batch to entire curve
censor.cv_batch = cell(nGridTime,1);
for t=1:nGridTime
    censor.cv_batch{t} = ones( 1, censor.nGrid(t)) * randi(censor.cv_K);
end

% calculate how much data points are in individual batches
censor.cv_batch_counts = zeros(1,censor.cv_K);
for kk = 1:censor.cv_K
    for t=1:censor.nGridTime
        censor.cv_batch_counts(kk) = censor.cv_batch_counts(kk) + sum(censor.cv_batch{t} == kk);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             fit of the mean function for X

[bw_mu_opt,mu_est_smoother_opt] = fCVforMU_periodic(censor,onb);
% [bw_mu_opt,mu_est_smoother_opt] = fKCVforMU(censor,onb);


mu_est_ONB = onb.smoothingMatrix_gridSmoother2ONB * mu_est_smoother_opt;
mu_est_inSpace = onb.onbMatrix*mu_est_ONB;

% visualize
if show_plots
    subplot(1,1,1)
    plot( mu_est_inSpace )
    pause(.01)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             fit covariance kernels


% disp('I am skipping CV for the surface smoothing')
bw_r = 0.1; % CV= 0.06336
bw_v = 0.2; % CV= 0.15672
% covLag0_smoother = fEstim_CovLagh(censor,onb,0,mu_est_inSpace,bw_r);

[covLag0_smoother, covLag0_ridge_smoother, sigma2_est_loc_qv, sigma2_est_loc_lin] = fnoCVforR_locQv(censor,onb,mu_est_inSpace,bw_r, bw_v);

covLag0_ONB = onb.smoothingMatrix_gridSmoother2ONB * covLag0_smoother * onb.smoothingMatrix_gridSmoother2ONB';
covLag0_inSpace = onb.onbMatrix * covLag0_ONB * onb.onbMatrix';

if show_plots
    subplot(1,1,1)
    surf( covLag0_inSpace )
    pause(.01)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       fit spectral density

nGridFreq = 1000;

%% estimate the spectral density
bw_F = bw_r; % the same smoothing parameter as for cov lag 0
specDensity_v7 =  fEstimateSpecDensity_v7_2020(censor,onb,mu_est_inSpace,bw_F,nGridFreq);
specDensity_v7s = fEstimateSpecDensity_v7_2020(censor,onb,mu_est_inSpace,0.15,nGridFreq);



%% extract auto-covariances (integrate out the spectral density)
numOfLags_all = censor.nGridTime - 1;
% f212_covLags = fEstimateCovLags(onb,f212_specDensity,numOfLags_all);
v7_covLags = fEstimateCovLags(onb,specDensity_v7,numOfLags_all);
v7s_covLags = fEstimateCovLags(onb,specDensity_v7s,numOfLags_all);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     kriging for the process X

% if I want to get confidence bands
cofidence_bands = 1;

% how much to padd at both ends - used for lagged regression
dataFull_padding = 30;

% use the sigma from the loc-qv estimation
s2 = sigma2_est_loc_qv;
if s2 < 1e-2, s2 = 1e-2; end
sigma_est = sqrt( s2 );

[krigingX_dyna,krigingX_dyna_padded] = fKrigingX( censor, onb, v7_covLags, mu_est_ONB, sigma_est, 1, cofidence_bands, dataFull_padding );
[krigingX_dyna_s,krigingX_dyna_padded_s] = fKrigingX( censor, onb, v7s_covLags, mu_est_ONB, sigma_est, 1, cofidence_bands, dataFull_padding );


save('Wank_matlab_variables_after_a.mat')
%% visualise kriging results


%% some plots

if show_plots
    nDays = 3*5;
    offset = 40 * 2;
    for ii = 1:nDays
        t = offset+ii;
        subplot(3,5,ii)
        plot( Eat_daily(t,:) )    
        hold on
        nonzero_positions = find(Eat_indic(t,:));
        scatter( nonzero_positions, Eat_daily(t, nonzero_positions), 'xb' )
%         plot( 100*NOx_daily(t,:) )
        hold off
        title("day="+num2str(t))
    end 

end

















