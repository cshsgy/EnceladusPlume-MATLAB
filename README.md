#EnceladusPlume-MATLAB

To run the code, MATLAB version of 2018b or higher may be required.

2022_Jun_23_Single_Func_Solver/composed_main_func.m is the function to run for the results shown
in the paper. The function takes three inputs: wmin, wmax, and depth, and will save the
results to a .mat data file "results.mat".

It runs quite slow right now (to prevent oscillation caused by an unstable wall surface
temperature), we are working on a C++ solution and will released it in the future.

Updated in 2023: Created a revised solver with easier access to most functions.
