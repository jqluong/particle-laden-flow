# Bidisperse Particle Laden Flow

This code supports the work "The 'Brazil-nut effect' in bidisperse particle laden flow on an incline" by Jack Luong, Sarah Burnett, and Andrea Bertozzi.
This code originated from Jeffery Wong supporting his paper ["A conservation law model for bidensity suspensions on an incline"](https://www.sciencedirect.com/science/article/pii/S0167278916302123) appearing in *Physica D: Nonlinear Phenomena*.

## Folder Structure

TODO: add specific files that are "main"

**bidisperse_pde**: Contains the code that solves the system of PDEs that model bidisperse particle laden flow on an incline at a given angle $\alpha$ where a corresponding flux table (.mat file) for $\alpha$ is present in the working directory. Such a flux table can be generated with the files in the "ODEsolve" folder.

**fluxtable_\***: Each folder refers to a specific set of particle pairs.  In each folder, the flux table also corresponds to an inclination angle $\alpha$. These fluxtables are required to simulate the system of PDEs modelling bidisperse particle laden flow. The fluxtables used for figures in this paper are in the folder "fluxtable_paper".

**ODEsolve**: Contains the code that solves the equilibirum profile boundary value problem. The function "set_constants.m" should be edited to reflect different physical parameters (such as liquid viscosity, particle diameter. etc.). The script "bidensity_ODE_tester.m" will solve the eqilibrium profile ODE system and plot its solution.
