Functional Lagged Regression with Sparse Noisy Observations
===========================================================

Supporting material for the paper:
	**Rubín, Tomáš, and Victor M. Panaretos.** *"Functional lagged regression with sparse noisy observations."* Journal of Time Series Analysis 41.6 (2020): 858-882.
	
	https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12551.


A. CONTENTS
-----------

This repository represents the code used for the methodology presented by Rubin and Panaretos (2019) for the treatment of the functional lagged regression problem where the input regressor time series is a sparsely and irregularly observed functional time series. The code is organised in the following way

- **master**
This folder includes all the functions used for the simulation, estimation, prediction, and regression treatments of sparsely observed functional time series. The functions saved in this folder are called from the other two folders as well as from the demo file.

- **simulation_study**
The code in this folder was used as the foundation for the simulation study presented in Rubin and Panaretos (2019).

- **Wank_analysis**
The code in this folder was used for the analysis of the meteorological data at the Wank mountain as presented in Rubin and Panaretos (2019).

- **demo.m**
This script presents a demonstration on how to use the repository. It simulates a sparsely observed functional time series as well as a response time series. The script estimates the spectral density, the cross-spectral density, and the filter coefficients, and performs the prediction of the response time series. All estimates are plotted by the MATLAB plot functions. Refer to the comments in this script for details on individual steps.
	
	

B. REQUIREMENTS
---------------

The scripts are written in MATLAB and run/testen in its version R2018a.
The simulation requires the package 'fdaM' (https://www.psych.mcgill.ca/misc/fda/). The packages is included in the folder for each simulation.

Note that the simulations were run on a computational cluster.

C. USAGE
--------

Individuals are free to use the codes for the purpose academic research, provided it is properly acknowledged. For any other use, permission must first be arranged with the author(s). Unless otherwise specified, the author of the codes is Tomas Rubin (tomas.rubin@gmail.com). Please contact me if you find errors in the codes.


D. CONTACTS
------------------
Tomas Rubin

tomas.rubin@gmail.com

https://www.linkedin.com/in/tomas-rubin/

www.tomasrubin.com

https://github.com/tomasrubin


E. REFERENCES
----------------

Principal reference:
	
	**Rubín, Tomáš, and Victor M. Panaretos.** *"Functional lagged regression with sparse noisy observations."* Journal of Time Series Analysis 41.6 (2020): 858-882.
	
	https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12551
			

How to cite the code:

	DOI 10.5281/zenodo.3190740.

	https://github.com/tomasrubin/sparse-functional-lagged-regression
	
	Tomas Rubin and Victor M. Panaretos. Github repository: Functional lagged regression with sparse noisy observations. https://doi.org/10.5281/zenodo.3190740, 2019.

