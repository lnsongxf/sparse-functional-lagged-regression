function launcher_1_full(job_id)

% initialized the rng based on the job id
rng(job_id);

% run once each of the settings
for simulated_process = [1 4]
    for response_process = [1 2 3]
        for nGridTime = [300 600 900]
            initial_seed = randi(2^32)-1;
            sim_1_full(initial_seed,simulated_process,response_process,nGridTime)
        end
    end
end

end