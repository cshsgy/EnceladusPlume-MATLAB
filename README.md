#EnceladusPlume-MATLAB

To run the code, MATLAB version of 2018b or higher may be required.

2022_Jun_23_Single_Func_Solver/composed_main_func.m is the function to run for the results shown
in the paper. The function takes three inputs: wmin, wmax, and depth, and will save the
results to a .mat data file "results.mat".

It runs quite slow right now, since MATLAB's for-loop is VERY slow, and we are taking very small time
steps (since the characteristic timescale for the first several wall grids is quite small), we are 
working on a C++ solution and will released it in the future.

Updated in 2023: Created a revised solver with easier access to most functions. We have two separate modules
for the convenience of individual use.
In Gas_Dynamics_Interpolator, the gas dynamics is pre-calculated with width and depth of the water level as the
input, and each time the two variables are supplied for a solution. This is way faster than solving the wall every
step, which is good for an estimation of the mass flux.
In Liquid_Dynamics_Solver, the liquid dynamics is studied with a PDE solver that solves the coupled equation of
velocity and the water level in response to the wall movement. 4-th order Runge-Kutta method and adaptive time steps
are used, in which an abrupt expansion/contraction can cause the solution to slow down significantly.
