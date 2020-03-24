Visibility data at Wank Mountain
--------------------------------

The code analyses the visibility data at Wank mountain. The code uses the packages saved in the "master" folder. The scripts should be run in the following order:

a_Wank_spectral_X.m
	This script runs the spectral analysis of the sparsely observed functional time series of the fair-weather atmospheric electricity. It estimates the spectral density and the complete second order covariance structure of the data. Furthermore the script performs the functional data recovery procedure to predict the latent functional data. All the MATLAB variables are saved in "Wank_vis_estimated_dynamics_X". Running this script takes long and requires a lot of memory. 

b_Wank_regression_on_E.m
	This script runs the lagged regression method based on sparsely observed functional time series of fair-weather atmospheric electricity. It estimates the filter coefficients by both Tikhonov and truncation regularization technique and predicts response on the test partition.

c_Wank_regression_on_T.m
	This script runs the lagged regression method based on fully functional time series of temperature data. It estimates the filter coefficients by both Tikhonov and truncation regularization technique and predicts response on the test partition.

d_Wank_regression_on_ET.m
	This script runs the lagged regression method based on both time series simultaneously, estimates the filter coefficients by both Tikhonov and truncation regularization technique and predicts response on the test partition. It produces the plots of the estimated filter coefficients used in the paper and saves them in the folder "figures_paper". It also writes into the console the table of prediction MSE and R2 coefficients of determination.

e_Wank_regression_on_ET_visualisation.m
	This script creates the building blocks for the schematic figure on how the lagged regression works used in the paper and saves it in the folder "figures_paper".


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