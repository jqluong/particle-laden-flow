% Compares simulations to experimental data.
%
% May need to be reconfigured if the format of the exp. data changes.
% 
%
% The solver is called as:
%
% bsol = bidisperse_PDE_exp(ang, phiL, XL, simulation_time, tdom, V_slurry, 
%                              x_trans, run_name, override)
%
% ang is in deg, V_slurry in ml, simulation_time in minutes and x_trans = [L1, L2] is in cm.
% tdom is the set of times at which the profile is to be saved (note:
% front positions are saved at each time step regardless).
%
% Set override to 1 to re-compute data, otherwise it tries to use the
% existing data file.(saved to dir/runs_autosave)
%
% The plotter is called as:
%
% plot_soln(bsol,comparison_name, mono)
%
% where mono = 0 for bidisperse, 1 for mono. with only lighter particles
% and 2 for mono with only heavier particles (this is important for
% consistency in plot styles).
%
% The plotter can be replaced with whatever works.

%TODO: ADD IN CORRECT VOLUMES

figure(1)
figure(2) %These need to exist for the solver to output plots per step.

tf = 20; %time in minutes to run the simulation
nout = 41; %number of output frames
t_out = linspace(0,tf,nout);

%------------

X = 0.5;
phis = [0.2 0.25 0.3];
angs = [15 20 25];

%                  ---angle--->
first_transient = [25 32 36;   %  |
                   27 30 33;   % phi
                   24 28 32 ]; %  V
               
second_transient = [20 20 20;
                    20 20 20;
                    20 20 20 ];
                
second_transient = first_transient; %disable second transient

actual_volume = [100 100 100;
                 100 100 100;
                 100 100 100 ];
             
sols = cell(3,3);

for i=1:3
    for j=1:3
        phi = phis(i); ang = angs(j);
        L1 = first_transient(i,j);
        L2 = second_transient(i,j);
        vol = actual_volume(i,j);
        
        vol = 80;
        
        rname = sprintf('50GSB5_%dvf_%ddeg',100*phi,ang);
        sols{i,j} = bidisperse_PDE_exp(ang, phi, X, tf, t_out, vol, [L1 L2], [rname,'_num'],0);
        
        plot_soln(sols{i,j},[rname,'_100mL'],0);
    end
end



