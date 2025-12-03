# Euro OIS Term Structure Sensitivity to ECB Policy Rate Surprises

This repository contains the MATLAB code used for the analysis presented in the paper:

Herzel, S., and Nicolosi, M., “Sensitivity of the Euro OIS Term Structure to ECB Policy Rate Surprises”, .....

If you use this code, please cite the original paper.

The purpose of this code is to calibrate the quantitative finance model described in the paper and to perform the sensitivity analysis of the Euro OIS term structure to ECB monetary policy surprises.

The main scripts included in this repository are:

1. script_calibration_QF.m

This script calibrates the model using the dataset stored in data_QF.mat.

The calibrated model parameters and results are saved in model_QF.mat.

2. script_analisiMP4.m

This script implements the monetary policy sensitivity analysis of the Euro OIS term structure, following the methodology described in the paper. It uses the calibrated model output stored in model_QF.mat as well as teh original data in "data_QF.mat".

Disclaimer

This code is provided for research and educational purposes only. No warranty is given regarding accuracy or completeness.
