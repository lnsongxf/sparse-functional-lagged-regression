%{

Run this file to produces graphs used in the paper for "forecasting error".
- it requires simulation results saved in the folder "results/"
- the figures are saved in the folder "figures_paper/"


Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fileName_Tikhonov = 'results/sim_tikhonov.csv';
fileName_truncation = 'results/sim_threshold_value.csv';
fileName_complete_Tikhonov = 'results/sim_fullyfunctional_tikhonov.csv';
fileName_complete_truncation = 'results/sim_fullyfunctional_threshold_value.csv';

% if ~exist('dataFull_truncation')
    dataFull_complete_truncation = readtable(fileName_complete_truncation);
    dataFull_complete_Tikhonov = readtable(fileName_complete_Tikhonov);
    dataFull_truncation = readtable(fileName_truncation);
    dataFull_Tikhonov = readtable(fileName_Tikhonov);
% end

[simulation_cases,simulation_runs_iCase,simulation_runs_process,simulation_runs_response_process,simulation_runs_response_process_shape] = fPrepareSimulationCases();


all_truncation = cell(2,1);
all_Tikhonov = cell(2,1);
all_oracle = cell(2,1);

response_process_shapes = [2 3];

for response_process_shape_id = 1:2
    response_process_shape = response_process_shapes(response_process_shape_id);

    data_truncation = dataFull_truncation( dataFull_truncation.response_process_shape == response_process_shape, : );
    data_Tikhonov = dataFull_Tikhonov( dataFull_Tikhonov.response_process_shape == response_process_shape, : );
    data_complete_truncation = dataFull_complete_truncation( dataFull_complete_truncation.response_process_shape == response_process_shape, : );
    data_complete_Tikhonov = dataFull_complete_Tikhonov( dataFull_complete_Tikhonov.response_process_shape == response_process_shape, : );

    process_name = ["iid", "FAR(1, 0.7)", "FMA(4)"];
    process_simulated = [0 1 4];

    %% prepare the results for fully functional cases

    % vector_fullyfunctional_Tikhonov_MSEs = zeros(3,1);
    % vector_fullyfunctional_truncation_MSEs = zeros(3,1);


    trunctation_MSE_per_case_fullyfunctional = nan(1,4);
    trunctation_MSE_per_case_fullyfunctional_n = zeros(1,4);
    Tikhonov_MSE_per_case_fullyfunctional = nan(1,4);
    Tikhonov_MSE_per_case_fullyfunctional_n = zeros(1,4);
    oracle_MSE_per_case_fullyfunctional = nan(1,4);
    
    
    for ii = 1:4
        T = 300 * ii;
        
        % truncation
        MSEs = data_complete_truncation.ZErrorTest( ...
            (data_complete_truncation.T == T) & ...
            (data_complete_truncation.simulated_process >= 1 ));
        
        if ~isempty(MSEs)
            trunctation_MSE_per_case_fullyfunctional(ii) = median(real(MSEs));
            trunctation_MSE_per_case_fullyfunctional_n(ii) = length(MSEs);
        end
        
        % oracle
        MSEs = data_complete_truncation.ZErrorTest_trueDynamics_( ...
            (data_complete_truncation.T == T) & ...
            (data_complete_truncation.simulated_process >= 1 ));
        
        if ~isempty(MSEs)
            oracle_MSE_per_case_fullyfunctional(ii) = median(real(MSEs));
        end
        
        % Tikhonov
        
        MSEs = data_complete_Tikhonov.ZErrorTest( ...
            (data_complete_Tikhonov.T == T) & ...
            (data_complete_Tikhonov.simulated_process >= 1));
        
        if ~isempty(MSEs)
            Tikhonov_MSE_per_case_fullyfunctional(ii) = median(real(MSEs));
            Tikhonov_MSE_per_case_fullyfunctional_n(ii) = length(MSEs);
        end
        
    end

    
    
    %% sparsely observed
    nCases = 4*4;
    
    truncation_MSE_per_case = nan(1,nCases);
    truncation_MSE_per_case_n = zeros(1,nCases);
    Tikhonov_MSE_per_case = nan(1,nCases);
    Tikhonov_MSE_per_case_n = zeros(1,nCases);
    oracle_MSE_per_case = nan(1,nCases);
    
    for iCase = 1:nCases
        
        % truncation
        MSEs = data_truncation.ZErrorTest( ...
            (data_truncation.iCase == iCase) & ...
            (data_truncation.simulated_process >= 1));
        
        if ~isempty(MSEs)
            truncation_MSE_per_case(iCase) = median(real(MSEs));
            truncation_MSE_per_case_n(iCase) = length(MSEs);
        end
        
        % oracle
        MSEs = data_truncation.ZErrorTest_trueDynamics_( ...
            (data_truncation.iCase == iCase) & ...
            (data_truncation.simulated_process >= 1));
        
        if ~isempty(MSEs)
            oracle_MSE_per_case(iCase) = median(real(MSEs));
        end
        
        % Tikhonov        
        MSEs = data_Tikhonov.ZErrorTest( ...
            (data_Tikhonov.iCase == iCase) & ...
            (data_Tikhonov.simulated_process >= 1));
        
        if ~isempty(MSEs)
            Tikhonov_MSE_per_case(iCase) = median(real(MSEs));
            Tikhonov_MSE_per_case_n(iCase) = length(MSEs);
        end        
    end


    all_truncation{response_process_shape_id} = [reshape(truncation_MSE_per_case,4,4), trunctation_MSE_per_case_fullyfunctional']';
    all_Tikhonov{response_process_shape_id} = [reshape(Tikhonov_MSE_per_case,4,4),Tikhonov_MSE_per_case_fullyfunctional']';
    all_oracle{response_process_shape_id} = [reshape(oracle_MSE_per_case,4,4),oracle_MSE_per_case_fullyfunctional']';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create graphs - dependence on the sample size
h_fig = figure('Position', [200,200, 1200, 600]);

% this is a nice graph
response_process_shape_names = {'A','B'};
T_names = {'300','600','900','1200'};
for shape_id = 1:2
    for T_id = 1:4
        subplot(2,4,T_id + (shape_id-1)*4 );
        plot( all_truncation{shape_id}(:,T_id), 'k-x' )
        hold on
        plot( all_Tikhonov{shape_id}(:,T_id), 'k:o' )
        plot( all_oracle{shape_id}(:,T_id), 'k--*' )
        hold off
        xticks(1:5)
        xticklabels({'10','20','40','60','inf'})
        xlim([ 0.6 5.4])        
        title("T="+T_names{T_id}+", shape="+response_process_shape_names{shape_id} )
        ylabel('RMSE, forecast')
        xlabel('N^{max}')
        
        if shape_id == 1 && T_id <= 2
            lgd = legend({'trunc','Tikh','oracle'},'Location', 'southwest');
            lgd.NumColumns = 1;
        else
            lgd = legend({'trunc','Tikh','oracle'},'Location', 'northeast');
            lgd.NumColumns = 1;
        end
        
        
        
        if shape_id == 1
            ylim([-0.1 1.1])
        else
            ylim([-0.01 0.3])
        end
        
    end
end

set(h_fig,'Units','Inches');
pos = get(h_fig,'Position');
set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_fig, 'figures_paper/forecast_rmse.pdf','-dpdf','-r0')
close(h_fig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dependence on the simulated process and the response process

%% simulated process

truncation_per_process = cell(2,1);
Tikhonov_per_process = cell(2,1);
oracle_per_process = cell(2,1);
truncation_per_response = cell(2,1);
Tikhonov_per_response = cell(2,1);
oracle_per_response = cell(2,1);

for response_process_shape_id = 1:2
    response_process_shape = response_process_shapes(response_process_shape_id);

    data_truncation = dataFull_truncation( dataFull_truncation.response_process_shape == response_process_shape, : );
    data_Tikhonov = dataFull_Tikhonov( dataFull_Tikhonov.response_process_shape == response_process_shape, : );
%     data_complete_truncation = dataFull_complete_truncation( dataFull_complete_truncation.response_process_shape == response_process_shape, : );
%     data_complete_Tikhonov = dataFull_complete_Tikhonov( dataFull_complete_Tikhonov.response_process_shape == response_process_shape, : );

    
    % per simulated process
    simulated_processes = [1 4];
    truncation_per_process{response_process_shape_id} = nan(length(simulated_processes),1);
    Tikhonov_per_process{response_process_shape_id} = nan(length(simulated_processes),1);
    oracle_per_process{response_process_shape_id} = nan(length(simulated_processes),1);
    for simulated_process_id = 1:2
        simulated_process = simulated_processes(simulated_process_id);
        
        % truncation    
        MSEs = data_truncation.ZErrorTest(data_truncation.simulated_process == simulated_process);
        truncation_per_process{response_process_shape_id}(simulated_process_id) = median(MSEs);
        
        % oracle
        MSEs = data_truncation.ZErrorTest_trueDynamics_(data_truncation.simulated_process == simulated_process);
        oracle_per_process{response_process_shape_id}(simulated_process_id) = median(MSEs);
        
        % Tikhonov
        MSEs = data_Tikhonov.ZErrorTest(data_Tikhonov.simulated_process == simulated_process);
        Tikhonov_per_process{response_process_shape_id}(simulated_process_id) = median(MSEs);
    end    
    
    % per response process
    response_processes = [1 2 3];
    truncation_per_response{response_process_shape_id} = nan(length(response_processes),1);
    Tikhonov_per_response{response_process_shape_id} = nan(length(response_processes),1);
    oracle_per_response{response_process_shape_id} = nan(length(response_processes),1);
    for response_process_id = 1:3
        response_process = response_processes(response_process_id);
        
        % truncation
        MSEs = data_truncation.ZErrorTest(data_truncation.response_process == response_process);
        truncation_per_response{response_process_shape_id}(response_process_id) = median(MSEs);
        
        % oracle
        MSEs = data_truncation.ZErrorTest_trueDynamics_(data_truncation.response_process == response_process);
        oracle_per_response{response_process_shape_id}(response_process_id) = median(MSEs);
        
        % Tikhonov
        MSEs = data_Tikhonov.ZErrorTest(data_Tikhonov.response_process == response_process);
        Tikhonov_per_response{response_process_shape_id}(response_process_id) = median(MSEs);
    end
    
end


h_fig = figure('Position', [200,200, 1200, 300]);

for response_process_shape_id = 1:2

    subplot(1,2, response_process_shape_id )
    plot( 1:2, truncation_per_process{response_process_shape_id} , 'k-x' )
    hold on
    plot( 1:2, Tikhonov_per_process{response_process_shape_id} , 'k:o')
    plot( 1:2, oracle_per_process{response_process_shape_id} , 'k--*')
    plot( 3:5, truncation_per_response{response_process_shape_id} , 'k-x' )
    plot( 3:5, Tikhonov_per_response{response_process_shape_id} , 'k:o')
    plot( 3:5, oracle_per_response{response_process_shape_id} , 'k--*')
    hold off
    xticks(1:5)
    xticklabels({'FAR(1)','FMA(4)','reg=1','reg=2','reg=3'})    
    xlim([0.6 5.4])
    ylabel('RMSE, forecast')
    title("shape="+response_process_shape_names{response_process_shape_id})        
    
    if response_process_shape_id == 1
        ylim([0 1])
        lgd = legend({'trunc','Tikh','oracle'},'Location','southwest');
    else
        ylim([0 0.2])
        lgd = legend({'trunc','Tikh','oracle'},'Location','southwest');
    end
    lgd.NumColumns = 3;
    
%     legend('boxoff')

end


set(h_fig,'Units','Inches');
pos = get(h_fig,'Position');
set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_fig, 'figures_paper/forecast_rmse_per_process.pdf','-dpdf','-r0')
close(h_fig)

