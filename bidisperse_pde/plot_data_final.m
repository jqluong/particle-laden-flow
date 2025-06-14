%% This is a script to only print a given frame of a run, assuming 
%  you have a sol struct loaded already.

%Setting up parameters needed to plot
n = length(sol.t_out);
V_slurry = sol.Vsl*10^6;
g = 9.8;
Vsl = (V_slurry)*(1e-6); %slurry volume in m^3;
Ltr = 1; %track length in m
angr = ang*pi/180; %angle in radians
wtr = (14)*1e-2; % track width in m
H = Vsl/(wtr*Ltr); %char. height H
t_dim = A.nul*Ltr/(H^2*g*sin(angr)); %time scaling (in s)

% Choose which frame you want here. If you want the final frame, it's
% easiest to put in "end" as the index instead of "i".
i = 2;

%Plot
plot(sol.X,sol.Y(:,:,end),'LineWidth',4);
xlabel('Track length (m)')
xlim([-0.1 1.3])
ylabel('$$h$$, $$h\phi_1$$, $$h\phi_2$$', 'interpreter', 'latex')
legend("Fluid", "Species 1 (Larger)", "Species 2 (Smaller)")
fontsize(40, "points")
colororder(["#000000" "#FF0000" "#20ADFF"])
%The order of the variables is fluid, larger, smaller
%colororder(["#000000" "#FF0000" "#0000FF"]) %black, red, blue
title(sprintf('t = %5f s',sol.t_out(i)*t_dim));
set(gca,'TickLength',[0.025 0.025],'LineWidth',2);