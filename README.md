# Trialhistory_lapses_EA 

## Overview
This repository contains the code necessary for generating figures related to the manuscript titled "Trial-history biases in evidence accumulation can give rise to apparent lapses."

The associated data for this manuscript can be found at dx.doi.org/10.6084/m9.figshare.24113793

## Usage Instructions

To reproduce the figures presented in the paper, follow these steps:

1. Clone this repository to your local machine using your preferred method.

2. Download the required data from the figshare link provided above.

3. Edit the following paths in the `master_startup.m` file:
   - Base directory (`path.base`)
   - Data directory (`path.data_save`)
   - Source data directory (`path.mat_save')
   - Directory for saving figures (`path.fig_save`)

4. Run the `make_figures.m` script to generate the main figures. Supplementary figure code can be found in their respective folders.

## Accumulator Model Fits

The fits to accumulator models were performed using code developed within the Julia package available at [https://github.com/Brody-Lab/PulseInputDDM](https://github.com/Brody-Lab/PulseInputDDM).

If you are interested in using this code or have any questions, please feel free to get in touch with us at diksha.gupta@ucl.ac.uk.

For any inquiries related to the manuscript or this repository, don't hesitate to reach out. 
