# HCVIntracellularModelling
HCVIntra.f90 is a Fortran code that is used for the stochastic simulations of the HCV intracellular infection, corresponding to the paper submitted to Scientific Reports. The code is written in Fortran 90 and can be compiled using gfortran. The -O4 flag should be used to optimise the executable output. This code as currently configured the treatment-free control case. The results are written to the file "fort.9".

On a desktop computer (~2Ghz processor) the HCVIntra.f90 code should take around 10 mins to run and uses <1 GB of RAM.

Solve_ODE_Model.m is a matlab code that uses ode45 to solve the ODE system HCVIntraModel.m which is the deterministic version of the stochastic model.

On a desktop computer (~2Ghz processor) the Solve_ODE_Model.m code should take around 40 mins to run and uses around 64 GB of RAM.
