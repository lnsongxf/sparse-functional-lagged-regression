Functional Lagged Regression with Sparse Noisy Observations
Tomas Rubin, Victor M. Panaretos
===========================================================

A. CONTENTS
-----------

This folder contains the following 5 folders:
1. results
2. script - exp covariance, fully functional
3. script - exp covariance, sparse
4. script - nrq covariance, fully functional
5. script - nrq covariance, sparse

The folders 2.-5. contain script for simulating the considered functional time series processes, applying the sparse observation protocol, and the subsequent estimation and prediction methodology. Each of these 4 folders contains a script 'launcher_x.m' where x=1,2,3,4. Running this script requires a random seed, for which we use the job_id number assign by the computational cluster. The script 'launcher_x.m' calls one simulation run for each of the considered simulation settings (sample size parameters N^max and T, the dynamics of the simulated processes, the lagged regression design). The script creates CSV files with the results of each simulation run (see below).

The 1. folder 'results' contains 8 csv files that have been transfered from the folders 2.-5. after running each laucher file 200 times (on a computational cluster):
- exp_sim_complete_threshold.csv
- exp_sim_complete_tikhonov.csv
- exp_sim_threshold_value.csv
- exp_sim_tikhonov.csv
- nrq_sim_complete_threshold.csv
- nrq_sim_complete_tikhonov.csv
- nrq_sim_tikhonov.csv
- nrq_sim_threshold_value.csv
The files:
- exp_analyse_csv_tables.m
- nrq_analyse_csv_tables.m

involve a script to go through the results in the CSV files and calculate the mean prediction/estimation errors. The scripts call 'prepare_simulation_cases.m' where the mapping from the 'iCase' variable to the sample size variables (N^max and T) is defined. Each of the two scripts 'exp_analyse_csv_tables.m', 'nrq_analyse_csv_tables.m' outputs matrices of outputs to be parse into 'LaTeX tables.xlsx' which prepares the LaTeX code used in the paper.

B. REQUIREMENTS
---------------

The scripts are written in MATLAB and run/testen in its version R2018a.
The simulation requires the package 'fdaM' (https://www.psych.mcgill.ca/misc/fda/). The packages is included in the folder for each simulation.

Note that the simulations were run on a computational cluster.



--------------------------------------------------------------------------------------------------------------------------------------

Tomas Rubin
EPFL SB MATHAA SMAT
MA B1 493 (Batiment MA)
Station 8
CH-1015 Lausanne
Switzerland

Email: tomas.rubin@epfl.ch

