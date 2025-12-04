# Euro OIS Term Structure Sensitivity to ECB Policy Rate Surprises

This repository contains the MATLAB code used for the analysis presented in the paper:

Herzel, S., and Nicolosi, M., “Sensitivity of the Euro OIS Term Structure to ECB Policy Rate Surprises”, ....

If you use this code, please cite the original paper.

The purpose of this code is to calibrate the model described in the paper and to perform the sensitivity analysis of the Euro OIS term structure to ECB monetary policy surprises.

The main scripts included in this repository are:

1. script_calibration.m

This script calibrates the model using the time-series of yields on EuroSTR stored in data.mat. 

The code allows to calibrate the model selecting a different tome window (rows 71-72) as well as a asubset of times to maturity (rows 81-82) 

The calibrated model parameters and results are saved in calibrated_model.mat.

2. script_sensitivity.m

This script implements the monetary policy sensitivity analysis of the Euro OIS term structure, following the methodology described in the paper. 

It uses the output from script_calibration stored in calibrated_model.mat as well as the original data in "data_QF.mat".

Disclaimer

This code is provided for research and educational purposes only. No warranty is given regarding accuracy or completeness.
