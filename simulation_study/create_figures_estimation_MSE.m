%{

Run this file to produces graphs used in the paper for "estimation error".
- it requires simulation results saved in the folder "results/"
- the figures are saved in the folder "figures_paper/"


Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}


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

FAR_truncation= cell(2,1);
FMA_truncation = cell(2,1);
all_truncation = cell(2,1);
FAR_Tikhonov = cell(2,1);
FMA_Tikhonov = cell(2,1);
all_Tikhonov = cell(2,1);

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


    trunctation_MSE_per_case_fullyfunctional = nan(3,4);
    trunctation_MSE_per_case_fullyfunctional_n = zeros(3,4);
    Tikhonov_MSE_per_case_fullyfunctional = nan(3,4);
    Tikhonov_MSE_per_case_fullyfunctional_n = zeros(3,4);

    for ii_process = 1:4
        for ii = 1:4
            T = 300 * ii;

            % truncation
            if ii_process <= 3
                MSEs = data_complete_truncation.transfer_MSE( ...
                    (data_complete_truncation.T == T) & ...
                    (data_complete_truncation.simulated_process == process_simulated(ii_process)));
            else
                MSEs = data_complete_truncation.transfer_MSE( ...
                    (data_complete_truncation.T == T) & ...
                    (data_complete_truncation.simulated_process >= 1 ));
            end

            if ~isempty(MSEs)
                trunctation_MSE_per_case_fullyfunctional(ii_process,ii) = median(real(MSEs));
                trunctation_MSE_per_case_fullyfunctional_n(ii_process,ii) = length(MSEs);
            end

            % Tikhonov
            if ii_process <= 3
                MSEs = data_complete_Tikhonov.transfer_MSE( ...
                    (data_complete_Tikhonov.T == T) & ...
                    (data_complete_Tikhonov.simulated_process == process_simulated(ii_process)));
            else
                MSEs = data_complete_Tikhonov.transfer_MSE( ...
                    (data_complete_Tikhonov.T == T) & ...
                    (data_complete_Tikhonov.simulated_process >= 1));
            end
            

            if ~isempty(MSEs)
                Tikhonov_MSE_per_case_fullyfunctional(ii_process,ii) = median(real(MSEs));
                Tikhonov_MSE_per_case_fullyfunctional_n(ii_process,ii) = length(MSEs);
            end
        end
    end



    %% sparsely observed
    nCases = 4*4;

    truncation_MSE_per_case = nan(3,nCases);
    truncation_MSE_per_case_n = zeros(3,nCases);
    Tikhonov_MSE_per_case = nan(3,nCases);
    Tikhonov_MSE_per_case_n = zeros(3,nCases);

    for ii_process = 1:4
        for iCase = 1:nCases

            % truncation
            if ii_process <= 3
                MSEs = data_truncation.transfer_MSE( ...
                    (data_truncation.iCase == iCase) & ...
                    (data_truncation.simulated_process == process_simulated(ii_process)));
            else
                MSEs = data_truncation.transfer_MSE( ...
                    (data_truncation.iCase == iCase) & ...
                    (data_truncation.simulated_process >= 1));
            end
            

            if ~isempty(MSEs) > 0
                truncation_MSE_per_case(ii_process, iCase) = median(real(MSEs));
                truncation_MSE_per_case_n(ii_process, iCase) = length(MSEs);
            end

            % Tikhonov
            if ii_process <= 3
                MSEs = data_Tikhonov.transfer_MSE( ...
                    (data_Tikhonov.iCase == iCase) & ...
                    (data_Tikhonov.simulated_process == process_simulated(ii_process)));
            else
                MSEs = data_Tikhonov.transfer_MSE( ...
                    (data_Tikhonov.iCase == iCase) & ...
                    (data_Tikhonov.simulated_process >= 1));
            end

            if ~isempty(MSEs) > 0
                Tikhonov_MSE_per_case(ii_process, iCase) = median(real(MSEs));
                Tikhonov_MSE_per_case_n(ii_process, iCase) = length(MSEs);
            end        
        end
    end



   
    FAR_truncation{response_process_shape_id} = [reshape(truncation_MSE_per_case(2,:),4,4), trunctation_MSE_per_case_fullyfunctional(2,:)']';
    FMA_truncation{response_process_shape_id} = [reshape(truncation_MSE_per_case(3,:),4,4), trunctation_MSE_per_case_fullyfunctional(3,:)']';
    all_truncation{response_process_shape_id} = [reshape(truncation_MSE_per_case(4,:),4,4), trunctation_MSE_per_case_fullyfunctional(4,:)']';

    FAR_Tikhonov{response_process_shape_id} = [reshape(Tikhonov_MSE_per_case(2,:),4,4),Tikhonov_MSE_per_case_fullyfunctional(2,:)']';
    FMA_Tikhonov{response_process_shape_id} = [reshape(Tikhonov_MSE_per_case(3,:),4,4),Tikhonov_MSE_per_case_fullyfunctional(3,:)']';
    all_Tikhonov{response_process_shape_id} = [reshape(Tikhonov_MSE_per_case(4,:),4,4),Tikhonov_MSE_per_case_fullyfunctional(4,:)']';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create graphs - dependence on the sample size
h_fig = figure('Position', [200,200, 1200, 500]);

% this is a nice graph
response_process_shape_names = {'A','B'};
T_names = {'300','600','900','1200'};
for shape_id = 1:2
    for T_id = 1:4
        subplot(2,4,T_id + (shape_id-1)*4 );
        plot( all_truncation{shape_id}(:,T_id), 'k-x' )
        hold on
        plot( all_Tikhonov{shape_id}(:,T_id), 'k:o' )
        hold off
        xticks(1:5)
        xticklabels({'10','20','40','60','inf'})
        xlim([ 0.7 5.3])
        ylim([0 1.1])
        title("T="+T_names{T_id}+", shape="+response_process_shape_names{shape_id} )
        ylabel('MSE, filters')
        xlabel('N^{max}')
        
        lgd = legend({'trunc','Tikh'},'Location','southwest');
        lgd.NumColumns = 2;
    end
end

set(h_fig,'Units','Inches');
pos = get(h_fig,'Position');
set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_fig, 'figures_paper/filter_mse.pdf','-dpdf','-r0')
close(h_fig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dependence on the simulated process and the response process

%% simulated process

truncation_per_process = cell(2,1);
Tikhonov_per_process = cell(2,1);
truncation_per_response = cell(2,1);
Tikhonov_per_response = cell(2,1);

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
    for simulated_process_id = 1:2
        simulated_process = simulated_processes(simulated_process_id);
        
        % truncation    
        MSEs = data_truncation.transfer_MSE(data_truncation.simulated_process == simulated_process);
        truncation_per_process{response_process_shape_id}(simulated_process_id) = median(MSEs);
        
        % Tikhonov
        MSEs = data_Tikhonov.transfer_MSE(data_Tikhonov.simulated_process == simulated_process);
        Tikhonov_per_process{response_process_shape_id}(simulated_process_id) = median(MSEs);
    end    
    
    % per response process
    response_processes = [1 2 3];
    truncation_per_response{response_process_shape_id} = nan(length(response_processes),1);
    Tikhonov_per_response{response_process_shape_id} = nan(length(response_processes),1);
    for response_process_id = 1:3
        response_process = response_processes(response_process_id);
        
        % truncation
        MSEs = data_truncation.transfer_MSE(data_truncation.response_process == response_process);
        truncation_per_response{response_process_shape_id}(response_process_id) = median(MSEs);
        
        % Tikhonov
        MSEs = data_Tikhonov.transfer_MSE(data_Tikhonov.response_process == response_process);
        Tikhonov_per_response{response_process_shape_id}(response_process_id) = median(MSEs);
    end
    
end



% for response_process_shape_id = 1:2
% 
%     subplot(1,4,1 + (response_process_shape_id-1)*2 )
%     plot( truncation_per_process{response_process_shape_id} , 'k-x' )
%     hold on
%     plot( Tikhonov_per_process{response_process_shape_id} , 'k:o')
%     hold off
%     xticks(1:2)
%     xticklabels({'FAR(1)','FMA(4)'})
%     ylim([0 1.2])
%     xlim([0.6 2.4])
%     ylabel('MSE, filters')
%     title("shape="+response_process_shape_names{response_process_shape_id})
% 
%     subplot(1,4,2 + (response_process_shape_id-1)*2 )
%     plot( truncation_per_response{response_process_shape_id} , 'k-x' )
%     hold on
%     plot( Tikhonov_per_response{response_process_shape_id} , 'k:o')
%     hold off
%     xticks(1:3)
%     xticklabels({'reg=1','reg=2','reg=3'})
%     ylim([0 1.2])
%     xlim([0.6 3.4])
%     ylabel('MSE, filters')
%     title("shape="+response_process_shape_names{response_process_shape_id})
% 
% end



h_fig = figure('Position', [200,200, 1200, 300]);

for response_process_shape_id = 1:2

    subplot(1,2, response_process_shape_id )
    plot( 1:2, truncation_per_process{response_process_shape_id} , 'k-x' )
    hold on
    plot( 1:2, Tikhonov_per_process{response_process_shape_id} , 'k:o')
    plot( 3:5, truncation_per_response{response_process_shape_id} , 'k-x' )
    plot( 3:5, Tikhonov_per_response{response_process_shape_id} , 'k:o')
    hold off
    xticks(1:5)
    xticklabels({'FAR(1)','FMA(4)','reg=1','reg=2','reg=3'})
    ylim([0 1.1])
    xlim([0.6 5.4])
    ylabel('MSE, filters')
    title("shape="+response_process_shape_names{response_process_shape_id})
    lgd = legend({'trunc','Tikh'},'Location','southwest');
    lgd.NumColumns = 2; 
%     legend('boxoff')

end


set(h_fig,'Units','Inches');
pos = get(h_fig,'Position');
set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_fig, 'figures_paper/filter_mse_per_process.pdf','-dpdf','-r0')
close(h_fig)

