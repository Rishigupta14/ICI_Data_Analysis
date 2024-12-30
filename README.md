# ICI_Data_Analysis
Related code from the course in Data Analysis and Uncertainty attended during the M.Sc. in Computational Fluid Dynamics at Cranfield University. The tasks proposed in the assignment are solved individually by a function named "f_taski.m" where (i) corresponds to the task number. The current state of the repository produces the figures as they appear in the report. However, the numerical values employed in it may not appear directly in the "read_data" workspace as they may have taken during a debugging operation.

On the other hand, the data obtained from the literature is stored in text files in the folder "Literature". The file should be loaded together with the functions for the overall post-processing software to work in nominal conditions.

To conclude, the visualization of the flow using ParaView is proportioned via the "3d_visualisation.pvsm" file. Note that no ParaView plots are included in the report by the author's decision. The state file was modified using ParaView's version 5.11.



## Overview
This project involves the analysis of velocity field data using various techniques including visualization, MATLAB scripting, experimental data acquisition, gradient computation, Fast Fourier Transformations (FFT), and statistical analysis.

## Repository Structure
- `3d_visualisation.pvsm`: ParaView state file for visualizing the instantaneous velocity field.
- `spatiotemporal_sample`: Folder containing the samples to be processed using MATLAB scripts.
- `read_data.m`: MATLAB file with comments explaining the structure of the data.

## Visualisation
The `3d_visualisation.pvsm` file is used to visualize the instantaneous velocity field. It employs specific detection methods to identify coherent structures. The state file ensures the reproduction of every visualization used in the project report.

## MATLAB Scripts
MATLAB scripts are created to process the samples provided in the `spatiotemporal_sample` folder. These scripts compute the mean value and the standard deviation of each velocity component at every sampling point. The results are plotted to identify the flow features responsible for the standard deviation peaks.

## Experimental Data Acquisition
The project involves a discussion on how a similar dataset could be acquired experimentally. It also evaluates the extent to which experimental and Direct Numerical Simulation (DNS) data are expected to overlap in canonical flows.

## Gradient Computation
The project computes the wall-normal gradient of the mean streamwise velocity component. An analytical solution is obtained based on symbolic operations in MATLAB. A grid convergence analysis and Richardson extrapolation are carried out based on uniform refinement.

## Fast Fourier Transformations (FFT)
FFT is carried out on the wall-normal velocity signal at the x=0.06 location. The effects of different window functions and various filtering techniques, such as moving averaging, are investigated. The frequency spectra are plotted and the frequency of the most energetic motions is identified.

## Statistical Analysis
The project investigates whether the wall-normal velocity component has a normal distribution by utilizing skewness and kurtosis measures and visualizing the probability density function of the velocity. A statistical test is employed to evaluate whether the correlation between any of the two velocity components is statistically significant.



