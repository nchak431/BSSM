# BSSM
Bayesian Spike and Slab Mixture model 

This contains all the R codes used for implementing the BSSM model proposed in the paper 'Bayesian Scalar-on-Image Regression with Spatial Interactions for Modeling Alzheimer's Disease' by N. Chakraborty, Q. Long and S. Kundu.

Details on the key files:

func.R

Contains all the R functions related to the estimation of proposed BSSM model. 

data_preparation.R

Code for preparing the brain image data obtained from ADNI-1 study to be used in our neuroimaging analysis. The code decomposes each quadrant of a particular 2-D slice of the T1w-MRI scans using wavelet decomposition.

predict_AD.R

Contains R code for applying BSSM model to the neuroimaging dataset obtained from ADNI-1 study and evaluate its prediction perfomrance.

generate_true_coeff.R

R code for generating true 2-D functional regression coefficients with different shapes including round, triangle and square, to be used in the simulations.

simulation.R

R code for performimg simuation studies using the BSSM model.

