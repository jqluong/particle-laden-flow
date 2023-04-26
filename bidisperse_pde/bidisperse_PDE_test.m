clear all
ang = 30;
A = set_constants2(ang); %This is where you change the angle

phi0 = 0.5;
X0 = 0.7;

tf = 20;
nout = 41;
t_out = linspace(0,tf,nout);

vol = 100;
L = [20 20];

figure(1)
figure(2) %These need to exist for the solver to output plots per step.

sol = bidisperse_PDE_exp(ang,phi0,X0,tf,t_out,vol,L,"ftableBD_30",1)