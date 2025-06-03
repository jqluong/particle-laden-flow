# Bidisperse Particle Laden Flow

This code supports the work "The 'Brazil-nut effect' in bidisperse particle laden flow on an incline" by Jack Luong, Sarah Burnett, and Andrea Bertozzi.
This code originated from Jeffery Wong supporting his paper ["A conservation law model for bidensity suspensions on an incline"](https://www.sciencedirect.com/science/article/pii/S0167278916302123) appearing in *Physica D: Nonlinear Phenomena*.
**Note:** Becuase this code originated from Jeffery Wong's work on bidensity particle laden flow, a lot of the function names still refer to the "bidipserse" model. This is only in name, the code correctly models the bidisperse case.

## Folder Structure

**bidisperse_pde**: Contains the code that solves the system of PDEs that model bidisperse particle laden flow on an incline at a given angle $\alpha$ where a corresponding flux table (.mat file) for $\alpha$ is present in the working directory. The flux tables can be generated with the script "generate_flux_table.m". It uses the functions defined in the folder "ODEsolve". Once a flux table is generated, the script "bidisperse_PDE_test.m" will simulate the system of conservation laws governing particle laden flow. 

**fluxtable_\***: Flux tables for various particle pairs are stored in these folders. In each folder, the flux tables also correspond to an inclination angle $\alpha$. These fluxtables are required to simulate the system of PDEs modelling bidisperse particle laden flow. The fluxtables used for figures in this paper are in the folder "fluxtable_paper".

**ODEsolve**: Contains the code that solves the equilibirum profile boundary value problem. The function "set_constants.m" should be edited to reflect different physical parameters (such as liquid viscosity, particle diameter. etc.). The script "bidensity_ODE_tester.m" will solve the eqilibrium profile ODE system and plot its solution.

**runs_autosave**: As a backup, the code will save the last used flux table here. This feature isn't used as much, but is left intact.

## How to Use

**If you want to recreate the figures from the paper**: For the equilibrium profiles, run the script "bidensity_ODE_tester.m" in folder "ODEsolve". For the full PDE simulations, make sure the flux table for the desired inclination angle is in the directory "bidisperse_pde". Then, run the script "bidisperse_PDE_test.m".

**If you want to run your own simulations**:
  1. Change the values in "set_constants.m" in the folder "ODEsolve" to your parameters.
  2. If you want equilibirum profiles, run the script "bidensity_ODE_tester.m" in folder "ODEsolve".
  3. If you want PDE simulations, you need to first generate the flux table using the script "generate_flux_table.m" in folder "bidisperse_pde". This takes about an hour.
  4. Double check the flux table is in the working directory. Run the script "bidisperse_PDE_test.m" in folder "bidisperse_pde". This takes about 10 minutes.
