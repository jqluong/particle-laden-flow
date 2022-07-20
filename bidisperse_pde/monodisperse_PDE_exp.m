% Solver for the monodisperse system to compare with experiments
% Jeffrey Wong
% Updated July 2014

% TODO: Merge with the bidisperse one to avoid code duplication

%-----Adjustable parameters-------

function bsol = monodisperse_PDE_exp(A, phiL, X0, simulation_time, V_slurry, x_trans, run_name, override)

ang = A.ang;
if(~override && exist([run_name,'.mat']))
    fprintf('Using existing data.\n');
    load(['./runs_autosave/',run_name,'.mat']);
    return
end

if(X0==1)
    rhos = A.rhos1;
else
    rhos = A.rhos2;
end

%CONVENTION: X = 0 or X = 1
%V_slurry = 84; %slurry volume in ml
particle_radius = 0.2; %particle radius in mm
%simulation_time = 20; %time to run the simulation in minutes

%nout = 24+1; %number of output frames (evenly spaced, includes t=0)
dx = 0.0025; %grid spacing

fprintf('\n Solving system with: \n phi0 = %f,\n angle = %f degrees\n\n',phiL,ang);

%-----Build fluxes prior to solve-----
%A = set_constants2(ang);

clear bsol;
bsol.angle = ang;
bsol.A = A;
%[fpars, dfpars] = get_fluxdata(fname); %Parameter list for fluxes
%J = @(u) bJac(u,fpars{:},dfpars{:});
%F = @(u) bF(u,fpars{:});

%------Define relevant parameters (and scalings)-------

angr = ang*pi/180; %angle in radians
%scaled setup
Vsl = (V_slurry)*(1e-6); %slurry volume in m^3;
wtr = (14)*1e-2; % track width in m
res_length = (7)*1e-2; %reservoir length in m
Ltr = 1; %track length in m
a = particle_radius/1e3; %particle radius in m
H = Vsl/(wtr*Ltr); %char. height H
g = 9.8; %gravity (m/s^2)

%------Compute some useful numbers------

%Crude estimate for transient time (seconds).
%t_trans = 9/2*A.nul*H/(g*A.rhos1*a^2*cos(angr)); 

%fprintf('First transient time: %f\n',t_trans_1);
%L_trans = 9/2*H^3*tan(angr)/(A.rhos1*a^2);

t_dim = A.nul*Ltr/(H^2*g*sin(angr)); %time scaling (in s)
fprintf('Time scale: %f\n',t_dim);

l_dim = Ltr; %length scale (in m)
u_dim = l_dim/t_dim; %velocity scale (in m/s)

%------More parameters
tf = simulation_time/(t_dim/60); %non. dim. time to run the simulation;

%-------Setup for ICs-------------

phiR = phiL; hR = 0.01; %small right state


%This computation isn't used (it's the huppert front estimate)

% C = (1 + A.rhos1*phiL*XL + A.rhos2*phiL*(1-XL));
% hL = Vsl/(res_length*wtr);
% C = C*(9*hL^2*res_length^2*9.8*sin(angr))/A.nul;
% C = C*(1-phiL/A.phimax)^2;
% C = C^(1/3);
% xf = C*t_trans_1^(1/3); %in m

res_length_eff = res_length;
rl = res_length_eff/Ltr;
hL = Ltr/res_length_eff; %non-dim

%fprintf('Effective reservoir length: %f\n',res_length_eff);

%----- Define the domain --------
xl = -1.1*rl; xr = 1.1*Ltr;
dt = dx*0.01; %initial timestep
%tdom = linspace(0,tf,nout);
X = xl:dx:xr; %domain
ngrid = numel(X);

%Processing initial conditions
yL = [hL hL*phiL];
yR = [hR hR*phiR]; 

%blob = @(x) +(-res_length/Ltr  <=x & x<= res_length_eff - res_length/Ltr);
blob = @(x) +(-rl <=x & x<= 0);
Y0 = ones(1,ngrid)'*yR + blob(X)'*(yL - yR);

%Save these into the solution struct for convenience
bsol.Y0 = Y0;
bsol.X = X;

%---------Run the actual solver-------------


set(0,'currentFigure',2);
plot(X,Y0(:,1),'-k',X,Y0(:,2),'--r');
title('Initial conditions');
legend('h','n');

CFL_max = 0.5;
Nx = length(X);
eps = 1e-12;
fronts = [1 1];
T = [0];
fluxes = zeros(Nx,2);
thresh = [hR + 0.02, phiR + 0.02]; %threshold for shock detection

Y = Y0; %remove this later

%density_factor = phiL*XL*A.rhos1 + phiL*(1-XL)*A.rhos2 + 1;
%nul_eff_1 = A.nul*(1-phiL*XL/A.phimax)^(-2)/density_factor;

% we solve 0 = h_t + C_wm*(h^3/3)_x and 0 = n_t + (C_wm h^3/3 n)_x,
% where the system is non-dim. by the scales for the full system.
% Here C_wm is the ratio of kinematic viscosities of the fluid to the
% fluid mixture.
C_wm = (phiL*rhos + 1)/(1-phiL/A.phimax)^(-2);

A.ODEsolve = @ode45;
[phi0s, fluxes_1, ~] = build_fluxtable1(A,X0,200);
mono_flux = @(Y) bF1_s(Y,phi0s,fluxes_1); %monodisperse flux function

t = 0;
it = 0;

dt_max = tf/100;
set(0,'currentFigure',1);

phase = 1;
max_speed = 0;
t_trans = 0; 
while(t < tf - eps)
    
    it = it + 1;
    T(it) = t;
    fprintf('%f and dt = %f, max speed = %f\n',t,dt,max_speed);
    
    for i=1:2
        itmp = find(Y(:,i) > thresh(i),1,'last');
        if(~isempty(itmp))
         fronts(it,i) = X(find(Y(:,i) > thresh(i),1,'last'));
        else
            fronts(it,i) = 0;
        end
    end
    
    if(phase==1 && (fronts(it,1) >= x_trans))
       fprintf('\n **********End of transient ************\n');
       t_trans = t*t_dim;
       phase = 2;
    end
    
    max_speed = 0;
    
    if(phase==1) %Well mixed (fluid + p1 + p2)
        for i=1:Nx
            fluxes(i,1) = C_wm*Y(i,1)^3/3; %WELL MIXED FLUXES
            fluxes(i,2) = C_wm*Y(i,2)*Y(i,1)^2/3;
            max_speed = max([max_speed,3*fluxes(i,:)/Y(i,1)]); %not quite a good estimate
        end
    else %p1 and p2 settled
        for i=1:Nx
            fluxes(i,:) = Y(i,1)^3*mono_flux(Y(i,:)');
            max_speed = max([max_speed,3*fluxes(i,:)/Y(i,1)]);
        end
    end
    
    %max_speed = 3*max(max(abs(fluxes./Y))); %estimate by max of f/h, g1/n1, g2/n2
    
    
    dt = CFL_max*dx/max_speed; %new dt
    dt = min(dt,dt_max);
    
    a = dt/dx;
    
    set(0,'currentFigure',1);
    
    plot(X,Y);
    title(sprintf('t = %f s',t*t_dim));
    drawnow
    
    %Y(1,:) = Y(1,:); %no change on left endpoint
    Ynew = Y(2:end,:) - a*(fluxes(2:end,:) - fluxes(1:(end-1),:));
    
    Y(2:end,:) = Ynew;
    %CFL estimate
    t = t + dt;
    %dt = 1e-4;
    if(t > tf - dt)
        dt = tf - t;
    end
end

fprintf('Transient time: %f\n',t_trans);


%Copy the rest of the results into the solution struct
bsol.T = T;
bsol.ang = ang;
bsol.Y = Y;
bsol.fronts = fronts;
%bsol.ictype = ic_type;
bsol.idata = [phiL hL; phiR hR];
bsol.Vsl = Vsl;
bsol.wtr = wtr;
bsol.res_length = res_length;
bsol.L = Ltr;
bsol.t_dim = t_dim;


%------------Save the data----------------

if(~isempty(run_name))
    save(['./runs_autosave/',run_name,'.mat'],'bsol');
end