# Bidisperse Particle Laden Flow

This code supports the work "The 'Brazil-nut effect' in bidisperse particle laden flow on an incline" by Jack Luong, Sarah Burnett, and Andrea Bertozzi.
This code originated from Jeffery Wong supporting his paper ["A conservation law model for bidensity suspensions on an incline"](https://www.sciencedirect.com/science/article/pii/S0167278916302123) appearing in *Physica D: Nonlinear Phenomena*.

## Folder Structure

TODO: add specific files that are "main"

**bidisperse_pde**: Contains the code that solves the system of PDEs that model bidisperse particle laden flow on an incline at a given angle $\alpha$ where a corresponding flux table (.mat file) for $\alpha$ is present in the working directory. The flux tables can be generated with the script "generate_flux_table.m". It uses the functions defined in the folder "ODEsolve". Once a flux table is generated, the script "bidisperse_PDE_test.m" will simulate the system of conservation laws governing particle laden flow. 

**fluxtable_\***: Flux tables for various particle pairs are stored in these folders. In each folder, the flux tables also correspond to an inclination angle $\alpha$. These fluxtables are required to simulate the system of PDEs modelling bidisperse particle laden flow. The fluxtables used for figures in this paper are in the folder "fluxtable_paper".

**ODEsolve**: Contains the code that solves the equilibirum profile boundary value problem. The function "set_constants.m" should be edited to reflect different physical parameters (such as liquid viscosity, particle diameter. etc.). The script "bidensity_ODE_tester.m" will solve the eqilibrium profile ODE system and plot its solution.
