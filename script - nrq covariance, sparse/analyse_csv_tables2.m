fileName_Tikhonov = 'z_sim_tikhonov.csv';
fileName_truncation = 'z_sim_threshold_value.csv';
fileName_complete_Tikhonov = 'z_sim_complete_tikhonov.csv';
fileName_complete_truncation = 'z_sim_complete_threshold_multiplier.csv';

dataFull_truncation = readtable(fileName_truncation);
dataFull_Tikhonov = readtable(fileName_Tikhonov);
dataFull_complete_truncation = readtable(fileName_complete_truncation);
dataFull_complete_Tikhonov = readtable(fileName_complete_Tikhonov);

data_truncation = dataFull_truncation;
data_Tikhonov = dataFull_Tikhonov;
data_complete_truncation = dataFull_complete_truncation;
data_complete_Tikhonov = dataFull_complete_Tikhonov;
% data_truncation = dataFull_truncation(dataFull_truncation.response_process == 1, :);
% data_Tikhonov = dataFull_Tikhonov(dataFull_Tikhonov.response_process == 1, :);


prepare_simulation_cases;

table_sim = nan(3,15+1+3,12);


% read files - sparse
simulated_processes = [0, 1, 4];
for simulated_process_id = 1:3
    simulated_process = simulated_processes(simulated_process_id);
    
    % tables by sample size
    for a = 1:3        
        for b = 1:5
            
            
            if b <= 4 % exclude the case with N^max = 40
                iCase = (b-1)*3 + a;
                
                % T first, then N
                row = 1 + (a-1)*5 + (b-1);
                table_sim(simulated_process_id, row, 1) = 300*a;
                table_sim(simulated_process_id, row, 2) = simulation_cases.nRandi_max(iCase);
                
                % first N, then T
                %             row = iCase;
                %             table_sim(simulated_process_id, row, 2) = 300*a;
                %             table_sim(simulated_process_id, row, 1) = simulation_cases.nRandi_max(iCase);
                
                
                subdata_truncation = data_truncation( (data_truncation.iCase == iCase) & (data_truncation.simulated_process == simulated_process), : );
                subdata_Tikhonov = data_Tikhonov( (data_Tikhonov.iCase == iCase) & (data_Tikhonov.simulated_process == simulated_process), : );
               

            else % complete observations
                row = 1 + (a-1)*5 + (b-1); 
                table_sim(simulated_process_id, row, 1) = 300*a;
                table_sim(simulated_process_id, row, 2) = 999;
                
                subdata_truncation = data_complete_truncation( (data_complete_truncation.T == 300*a) & (data_complete_truncation.simulated_process == simulated_process), : );
                subdata_Tikhonov = data_complete_Tikhonov( (data_complete_Tikhonov.T == 300*a) & (data_complete_Tikhonov.simulated_process == simulated_process), : );                               
                
            end            
            
            % transfer estimation
            table_sim(simulated_process_id, row,3) = mean( subdata_truncation.transfer_RMSE );
            table_sim(simulated_process_id, row,4) = var(  subdata_truncation.transfer_RMSE )^(1/2);
            table_sim(simulated_process_id, row,5) = mean( subdata_Tikhonov.transfer_RMSE );
            table_sim(simulated_process_id, row,6) = var(  subdata_Tikhonov.transfer_RMSE )^(1/2);
            
            % prediction
            table_sim(simulated_process_id, row,7) = mean( subdata_truncation.ZErrorTrain_holdout_ );
            table_sim(simulated_process_id, row,8) = var(  subdata_truncation.ZErrorTrain_holdout_ )^(1/2);
            table_sim(simulated_process_id, row,9) = mean( subdata_Tikhonov.ZErrorTrain_holdout_ );
            table_sim(simulated_process_id, row,10)= var(  subdata_Tikhonov.ZErrorTrain_holdout_ )^(1/2);
            table_sim(simulated_process_id, row,11)= mean( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
            table_sim(simulated_process_id, row,12)= var(  subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ )^(1/2);
        end
    end
    
    % tables by regression model
    for response_process = 1:3
        row = 16 + response_process;
        
        subdata_truncation = data_truncation( (data_truncation.response_process == response_process) & (data_truncation.simulated_process == simulated_process), : );
        subdata_Tikhonov = data_Tikhonov( (data_Tikhonov.response_process == response_process) & (data_Tikhonov.simulated_process == simulated_process), : );
        
        % transfer estimation
        table_sim(simulated_process_id, row,3) = mean( subdata_truncation.transfer_RMSE );
        table_sim(simulated_process_id, row,4) = var(  subdata_truncation.transfer_RMSE )^(1/2);
        table_sim(simulated_process_id, row,5) = mean( subdata_Tikhonov.transfer_RMSE );
        table_sim(simulated_process_id, row,6) = var(  subdata_Tikhonov.transfer_RMSE )^(1/2);
        
        % prediction
        table_sim(simulated_process_id, row,7) = mean( subdata_truncation.ZErrorTrain_holdout_ );
        table_sim(simulated_process_id, row,8) = var(  subdata_truncation.ZErrorTrain_holdout_ )^(1/2);
        table_sim(simulated_process_id, row,9) = mean( subdata_Tikhonov.ZErrorTrain_holdout_ );
        table_sim(simulated_process_id, row,10)= var(  subdata_Tikhonov.ZErrorTrain_holdout_ )^(1/2);
        table_sim(simulated_process_id, row,11)= mean( subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ );
        table_sim(simulated_process_id, row,12)= var(  subdata_Tikhonov.ZErrorTrain_holdoutTrueDynamics_ )^(1/2);
                
    end
    
end

squeeze(table_sim(1,:,:))
squeeze(table_sim(2,:,:))
squeeze(table_sim(3,:,:))
