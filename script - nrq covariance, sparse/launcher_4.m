function launcher_4(job_id)

% initialized the rng based on the job id
rng(job_id);

% run once each of the settings
for iCase = randperm(12)
    for simulated_process = [1 4]
        for response_process = [1 2 3]
            initial_seed = randi(2^32)-1;
            sim_91(initial_seed,simulated_process,response_process,iCase)
        end
    end
end

end