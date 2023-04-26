clear all
ang = 30;
A = set_constants2(ang); %This is where you change the angle

phi0 = 0.2;
X0 = 0.6;

tf = 20;
nout = 41;
t_out = linspace(0,tf,nout);

vol = 100; %Change this to 82.5
L = [0 0];

figure(1)
figure(2) %These need to exist for the solver to output plots per step.

sol = bidisperse_PDE_exp(ang,phi0,X0,tf,t_out,vol,L,"ftableBD_30",1)