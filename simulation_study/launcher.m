%{

This function is launches simulation runs. The individual simulation cases
have been split and grouped into "runs" given by "iRun" identifier so it
can be better fed into the computation cluster. "job_id" is the identifier
assigned by the cluster, which thanks to its uniqueness is used as random
seed. This script should be called by the BATCH command saved in the text
file "batch_commands_for_cluster/_sim_execute_jobs.txt"

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}


function launcher(job_id, iRun)

addpath('../master')
addpath('../master/fdaM')

% initialized the rng based on the job id
rng(job_id);

[simulation_cases,simulation_runs_iCase,simulation_runs_process,simulation_runs_response_process,simulation_runs_response_process_shape] =...
    fPrepareSimulationCases();

% run once each of the settings
for iCase = simulation_runs_iCase{iRun}(randperm(length( simulation_runs_iCase{iRun} )))
    for simulated_process = simulation_runs_process{iRun}(randperm(length(  simulation_runs_process{iRun} )))
        for response_process_shape = simulation_runs_response_process_shape{iRun}(randperm(length(  simulation_runs_response_process_shape{iRun} )))
            for response_process = simulation_runs_response_process{iRun}(randperm(length(  simulation_runs_response_process{iRun} )))
                initial_seed = randi(2^32)-1;
                simulation_run(initial_seed,simulated_process,response_process_shape,response_process,iCase)
            end
        end
    end
end


end