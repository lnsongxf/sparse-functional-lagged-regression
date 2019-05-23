

initial_seed=1
simulated_process=1
response_process=1
iCase = 4
prepare_simulation_cases;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the onb.basis

% nGridSpace = 241, nKnots = 24
onb = fPrepareONB(241, 24);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                          SIMULATION STARTS HERE                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% random seed that has been saved
tic


nGridTime = simulation_cases.nGridTime(iCase);
nRandi_max = simulation_cases.nRandi_max(iCase);


fprintf('#################################################################################\n')
fprintf('#################################################################################\n')
fprintf('\n')
fprintf('iCase = %d\n',iCase)
fprintf('seed = %d\n',initial_seed)
fprintf('nGridTime = %d\n',nGridTime)
fprintf('nRandi_max = %d\n',nRandi_max)
fprintf('simulated_process = %d\n',simulated_process)
fprintf('response_process = %d\n',response_process)
fprintf('\n')
fprintf('#################################################################################\n')
fprintf('#################################################################################\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate the process - twice
nGridFreq = 1000;
dataFull_padding = 10;

% if simulated_process <= 2
%     % AR(1)
%     pars_simulated_process = fSimulate_AR_random_prepare_coef(onb,simulated_process);
%     [censor, dataFull,dataFull_padded, covsLags_true, mu_true, specDensity_true, maxEigenvalue,sigma] = fSimulate_AR_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
%     [censor_test, dataFull_test,dataFull_test_padded,~,~,~,~,~] = fSimulate_AR_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
% else
%     % MA(q)
%     pars_simulated_process = fSimulate_MA_random_prepare_coef(onb);
%     [censor, dataFull,dataFull_padded, covsLags_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
%     [censor_test, dataFull_test,dataFull_test_padded,~,~,~,~,~] = fSimulate_MA_random(onb,simulated_process,pars_simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
% end


if simulated_process <= 2
    % AR(1)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, maxEigenvalue,sigma] = fSimulate_AR(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_AR(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
else
    % MA(q)
    [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covsLags_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
    [censor_test, dataFull_test,dataFullONB_test,dataFull_test_padded,dataFullONB_test_padded,~,~,~,~,~] = fSimulate_MA(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding);
end

%% set up the response process
pars_response_process = fSimulate_response_prepare_coef(onb.nGridSpace,10); % prepare 10 random functions
[response, transfer_true, specCrossDensity_true] = fSimulate_response_ONB(onb,response_process,pars_response_process,dataFullONB_padded,dataFull_padding,specDensity_true);
[response_test,~,~] = fSimulate_response_ONB(onb,response_process,pars_response_process,dataFullONB_test_padded,dataFull_padding,specDensity_true);

%% visualization

for t = 1:12
    subplot(3,4,t)
    plot(onb.gridSpace, dataFull(t,:))
    hold on
    plot( onb.gridSpace(censor.grid{t}) , censor.data{t},'x' )
    hold off
    integral = mean(dataFull(t,:));
    ylim([-4 4])
    xlim([0 1])
    title( ['t=',num2str(t),', int=',num2str(round(integral,2)),', Y=',num2str(round(response(t),2))] )
end
suptitle('first 12 snapshots')
pause(0.1)

xxx

%%

t=5
subplot(1,2,1)
plot(onb.gridSpace, dataFull(t,:), '-k')
ylim([-4 4])
xlabel('Domain x\in[0,1]')
ylabel('Value of X_t(x)')
title('Dense sampling')

subplot(1,2,2)
plot(onb.gridSpace, dataFull(t,:), ':k')
hold on
plot( onb.gridSpace(censor.grid{t}) , censor.data{t},'xk' )
hold off
ylim([-4 4])
xlabel('Domain x\in[0,1]')
ylabel('Value of X_t(x)')
title('Sparse sampling')

%%
subplot(2,1,1)
surf( real( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )
subplot(2,1,2)
surf( imag( onb.onbMatrix *transfer_true.specTransfer' ), 'EdgeColor','none' )

xxx
