%%This is a script used to solve the equilibrium system of ODEs.
% Requires: none
% Output: Figures that show the solutions to the equilibrium system.


clear all

%%Prompt user to input parameters: 

%Set angle
alpha = input("Enter in inclination angle of track (degrees): ");
A = set_constants2(alpha); %This is where you change the angle

%Prompt user for initial conditions
phi0 = input("Enter in particle fraction integral condition (phi_0): ");
X0 = input("Enter in particle ratio integral condition (chi_0): ");
p0 = [phi0 X0];  %Forms the initial guess for the shooting method. Using
%phi0 and X0 as the initial guess seems to work.

%%Solve ODE
sol = solve_bidensity_ODE(phi0,X0,p0,A);
Z = sol.Z;
X = sol.X;
phi = sol.phi;
sigma = sol.sigma;
T = sol.T;

%%Display Results
figure
plot(Z,phi,'-k','LineWidth', 2);
title("\phi v. Z for \phi_0: " + phi0 + " and \chi_0: " + X0 + " and \alpha: " + A.ang)
xlabel("z")
ylabel("\phi")

figure
plot(Z,X,'-k','LineWidth', 2);
title("\chi v. Z for \phi_0: " + phi0 + " and \chi_0: " + X0 + " and \alpha: " + A.ang)
xlabel("z")
ylabel("\chi")

figure
hold on
xlim([0,1])
plot(Z,X .* phi, 'LineWidth', 3, 'DisplayName', 'Larger Species')
plot(Z, (1-X) .* phi, 'LineWidth', 3, 'DisplayName', 'Smaller Species')
legend
title("\phi_1, \phi_2 v. z for \phi_0: " + phi0 + " and \chi_0: " + X0 + " and \alpha: " + A.ang)
xlabel("z")
ylabel("\phi")
fontsize(40, "points")
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
colororder(["#FF0000" "#20ADFF"]) %Variable order: larger, smaller
hold off


%%Debugging
%This was some debugging code used to check whether the integral conditions
%were being met

%{
trapz(Z,phi)
trapz(Z,X.*phi)/phi0
%}
