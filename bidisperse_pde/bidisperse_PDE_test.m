clear all
%Add a path
parentDirectory = fileparts(cd);
addpath(parentDirectory + "/ODEsolve")


%Set Constants
ang = 35;
A = set_constants2(ang); %This is where you change the angle

phi0 = 0.4;
X0 = 0.7;

tf = 20;
nout = 41;
t_out = linspace(0,tf,nout);
t_out = [t_out(1) t_out(2:end)];

vol = 82.5; %Change this to 82.5
L = [20 20];

figure(1)
figure(2) %These need to exist for the solver to output plots per step.

sol = bidisperse_PDE_exp(ang,phi0,X0,tf,t_out,vol,L,"ftableBD_35",1)