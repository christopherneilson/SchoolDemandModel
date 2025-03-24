# NFP Fake Data Testing Framework

This repository contains MATLAB scripts and functions for testing the NFP model described in Neilson (2025). The project focuses on simulating and evaluating non-linear market share models using fake data, optimized for debugging and validation.

## Setup Instructions

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/<your-username>/SchoolDemandModel.git

2. Ensure MATLAB is installed on your system. The repository is tested with MATLAB 2024b but should work with other recent versions.

3. Add the folder to MATLAB’s path:
   ```matlab
   addpath('/path/to/SchoolDemandModel');

#### **How to Run**

## How to Run

The `testing_NFP_simulatedData.m` script is the main testing file. It validates various components of the NFP model using simulated data. 

### Steps:
1. Open MATLAB and navigate to the repository folder.
2. Run the script:
   ```matlab
   testing_NFP_simulatedData

3. The script performs the following actions:
   - Loads (or creates) simulated data.
   - Checks the accuracy of the share equation.
   - Evaluates the inversion process for market shares.
   - Tests gradient computations.
   - Benchmarks optimization routines.

4.	Results are printed to the MATLAB command window, including parameter estimates, optimization performance, and diagnostic outputs.

## Workflow and Components

The script follows these steps:

### 1. Setting Parameters
- Global variables and parameters are initialized.

### 2. Data Handling
- Simulated data is loaded (`TestDataNFP_{date}.mat`). This is the last simulated dataset created.

### 3. Testing Modules
- **Share Equation Check**: Validates market share equations using `test_checkshareEquation.m`.
- **Inversion Process**: Evaluates the model's inversion capabilities using `test_inversion.m`.
- **Clocking Core Functions**: Benchmarks computational efficiency using `test_clockmainfunctions.m`.
- **Gradient Evaluation**: Tests the accuracy of gradient computations via `test_gradient_components.m`.

### 4. Optimization
- The script benchmarks optimization routines (`fmincon` or `knitromatlab`) for parameter estimation.
- Results include estimated coefficients, variance-covariance matrices, and objective function values.

### 5. Output
- Diagnostic tables (`mprint.m`) display parameter comparisons (True vs. Estimated).
- Performance metrics are saved for analysis.

## Supporting Files and Functions
- invertmarketshares.m: Computes market share inversion.
- gmmObj_NFP.m: Objective function for GMM estimation.
- getVarCov.m: Calculates variances and covariances.
- printMoments.m: Outputs model moments.
- test_checkshareEquation.m: Validates share equations.
- test_inversion.m: Tests the inversion process.
- test_clockmainfunctions.m: Benchmarks computational functions.
- test_gradient_components.m: Evaluates gradient accuracy.
- mprint.m: Prints formatted tables for parameter comparisons.
