%% set some constant concerning the simulations
nRandi_min = 0; % the minimal number of points per curve - the same for all cases
simulation_variables.nGridTime = [300 600 900]; % maximum time horizon T
simulation_variables.nRandi_max = [5 10 20 30 40]; % maximum sampled points per curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  prepare the structure for different simulation cases

nCases = length(simulation_variables.nGridTime) * length(simulation_variables.nRandi_max);

simulation_cases = [];
simulation_cases.nGridTime = zeros(1,nCases);
simulation_cases.nRandi_max = zeros(1,nCases);

ii = 1;
for i1 = simulation_variables.nRandi_max
    for i2 = simulation_variables.nGridTime    
            simulation_cases.nRandi_max(ii) = i1;
            simulation_cases.nGridTime(ii) = i2;            
            ii = ii + 1;
    end
end