Simulation study
--------------------------------

The code in this folder performs the simulation study presented in the paper. The simulation was run on a computational cluster and the individual runs were executed by BATCH commands saved in the file  "batch_commands_for_cluster/_sim_execute_jobs.txt". The results of the simulation runs are saved in the folder "results". These saved raw results are then accessed by "create_figures_estimation_MSE.m" and "create_figures_forecasting.m" to produce the figures for the paper.


create_figures_estimation_MSE.m
	This script generates the figures comparing the lagged regression filter coefficients estimates by MSE. The script outputs the figures into the folder "figures_paper".

create_figures_forecasting.m
	This script generates the figures comparing the prediction error. The script outputs the figures into the folder "figures_paper".

launcher.m
	This is a function handle that is called by a batch commend and groups together a few runs -> it calls the script simulation_run.m a few times with the designated setting.

prepare_batch_files.m
	This script creates the BATCH commands files saved in the folder "batch_commands_for_cluster".

simulation_run.m
	This script is the primary file handling the simulation. One runtime of this script runs a simulation of sparsely observed as well as fully observed functional time series with the designated setting and performs the estimation and prediction exercise. The results are saved into the CSV tables in the folder "results".


-----------------------------------
Developed and tested on: MATLAB R2018a


-----------------------------------
Individuals are free to use the codes for the purpose academic research, provided it is properly acknowledged. For any other use, permission must first be arranged with the author(s). Unless otherwise specified, the author of the codes is Tomas Rubin (tomas.rubin@gmail.com). Please contact me if you find errors in the codes.

-----------------------------------
Principal reference:
	
	Tomas Rubin and Victor Panaretos (2019). Functional Lagged Regression with Sparse Noisy Observations. arXiv:1905.07218.
	

The possibilities how to reference the code:

	DOI 10.5281/zenodo.3190740.

	github.com/tomasrubin/sparse-functional-lagged-regression
	
	Tomas Rubin and Victor M. Panaretos. Github repository: Functional lagged regression with sparse noisy observations. http://doi.org/10.5281/zenodo.3190741, 2019

	
	