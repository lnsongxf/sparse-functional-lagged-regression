%{

This script creates the BATCH files to be run by the computation cluster.
The command to be put into the console is saved in
"batch_commands_for_cluster/_sim_execute_jobs.txt".

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}



sim_id = '09'

[simulation_cases,simulation_runs_iCase,simulation_runs_process,simulation_runs_response_process,simulation_runs_response_process_shape] =...
    fPrepareSimulationCases();

nRuns = length(simulation_runs_iCase);
for iRun = 1:nRuns
    nGridTime = simulation_cases.nGridTime(max(simulation_runs_iCase{iRun}));
    nRandi_max = simulation_cases.nRandi_max(max(simulation_runs_iCase{iRun}));
    
    fileID = fopen(['batch_commands_for_cluster/la_',sim_id,'_run',num2str(iRun),'.txt'],'w');
    
    fprintf(fileID,'#!/bin/bash\n');

    fprintf(fileID,'#SBATCH --nodes 1\n');
    fprintf(fileID,'#SBATCH --ntasks 1\n');
    fprintf(fileID,'#SBATCH --time 5:59:00\n');
    fprintf(fileID,"#SBATCH --mem "+ round(simulation_cases.memory(simulation_runs_iCase{iRun}))*1024 +"\n");    
    fprintf(fileID,'#SBATCH --cpus-per-task=1\n');
    
    fprintf(fileID,['#SBATCH -o /home/trubin/sim_',sim_id,'/zmyjob.%%j.%%N.out\n']);
    fprintf(fileID,['#SBATCH -D /home/trubin/sim_',sim_id,'/\n']);

    fprintf(fileID,'module purge \n');
    fprintf(fileID,'module load matlab\n');
    fprintf(fileID,['matlab -nodisplay -r "launcher_',sim_id,'(${SLURM_JOB_ID}, ',num2str(iRun),')"\n']);
        
    fclose(fileID);
end


% prepare the BATCH command file
str = "for i in {1..90}; do ";
for iRun = flip(1:nRuns)
    for ii = 1:1    
        str = str + 'sbatch sim_'+sim_id+'/la_'+sim_id+'_run' + num2str(iRun) + '.txt; ';    
    end
end
str = str + "done";


fileID = fopen('batch_commands_for_cluster/_sim_execute_jobs.txt','w');
fprintf(fileID,str);
fclose(fileID);


