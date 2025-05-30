%%This is a script used to simulate the system of conservation laws that
% govern the evolution of bidisperse particle laden flow on a slurry.
% Requires: a flux table (some .mat file usually named ftable...)
% in the directory at the desired angle.
% Output: a .gif showing the time evolution of the system in this 
% directory.

clear all

%Add a path
parentDirectory = fileparts(cd);
addpath(parentDirectory + "/ODEsolve")

%%Set and prompt user for conditions
%Set angle
ang = 20; %Set the angle here. You should have a flux table at this angle.
A = set_constants2(ang);

%Prompt user for initial conditions
phi0 = input("Enter in initial total particle fraction (phi_0): ");
X0 = input("Enter in initial particle ratio (chi_0): ");

tf = 20;
nout = 41;
t_out = linspace(0,tf,nout);
t_out = [t_out(1) t_out(2:end)];

vol = input("Enter in initial total volume of mixture: ");
L_trans = input("Enter in transient length: "); %This is the length of the
%"well mixed" regime. The slurry will obey Huppert's 1/3 law and remain
%well mixed until it reaches L_trans;
L = [L_trans L_trans]; %The second element in the vector always
% should be the same as the first. We are not considering any "secondary"
% transient behavior in this work.

%%Simulate PDE
figure(1)
figure(2) %These need to exist for the solver to output plots per step.

%Get required flux table
ftablename = "ftableBD_" + num2str(ang);

%Call to solving the PDE
sol = bidisperse_PDE_exp(ang,phi0,X0,tf,t_out,vol,L,"ftablename",1);