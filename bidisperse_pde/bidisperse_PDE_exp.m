% Solver for the bidisperse system
% to compare with experiments
% Jeffrey Wong
% Updated July 2014

%-----Adjustable parameters-------

function bsol = bidisperse_PDE_exp(ang, phiL, XL, simulation_time, t_out, V_slurry, x_trans, run_name, override)

%ang = 15;
if(~override && exist([run_name,'.mat']))
    fprintf('Using existing data.\n');
    load([run_name,'.mat']);
    return
end
fname = ['ftable_',num2str(ang)];
%run_name = 'settled_15_ceramic_1'; %set this to a string to save the run
%phiL = 0.3; XL = 0; %reservoir concentration

x_trans_1 = x_trans(1)/100; %INPUT IS IN CM!
x_trans_2 = x_trans(2)/100; %set to negative number to skip the second transient

%V_slurry = 84; %slurry volume in ml
particle_radius = 0.2; %particle radius in mm
%simulation_time = 20; %time to run the simulation in minutes

%nout = 24+1; %number of output frames (evenly spaced, includes t=0)
dx = 0.0025; %grid spacing

fprintf('\n Solving system with: \n phi0 = %f, X0 = %f\n angle = %f degrees\n\n',phiL,XL,ang);

%-----Build fluxes prior to solve-----
A = set_constants2(ang);

clear bsol;
bsol.angle = ang;
bsol.A = A;
[fpars, dfpars] = get_fluxdata(fname); %Parameter list for fluxes
%J = @(u) bJac(u,fpars{:},dfpars{:});
F = @(u) bF(u,fpars{:});

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
t_out = t_out/t_dim*60; %input t_out is in minutes; convert it
fprintf('Time scale: %f\n',t_dim);

l_dim = Ltr; %length scale (in m)
u_dim = l_dim/t_dim; %velocity scale (in m/s)

%------More parameters
tf = simulation_time/(t_dim/60); %non. dim. time to run the simulation;

%-------Setup for ICs-------------

phiR = phiL; XR = XL; hR = 0.01; %small right state


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
xl = -1.2*rl; xr = 1.2*Ltr;
dt = dx*0.01; %initial timestep
%tdom = linspace(0,tf,nout);
X = xl:dx:xr; %domain
ngrid = numel(X);

%Processing initial conditions
yL = [hL hL*phiL*XL hL*phiL*(1-XL)];
yR = [hR hR*phiR*XR hR*phiR*(1-XR)]; 

%blob = @(x) +(-res_length/Ltr  <=x & x<= res_length_eff - res_length/Ltr);
blob = @(x) +(-rl <=x & x<= 0);
Y0 = ones(1,ngrid)'*yR + blob(X)'*(yL - yR);

%Save these into the solution struct for convenience
bsol.Y0 = Y0;
bsol.X = X;

%---------Run the actual solver-------------


set(0,'currentFigure',2);
plot(X,Y0(:,1),'-k',X,Y0(:,2),'--r',X,Y0(:,3),'--b');
title('Initial conditions');
legend('h','n1','n2');

CFL_max = 0.5;
Nx = length(X);
eps = 1e-12;
fronts = [0 0 0];
T = [0];
fluxes = zeros(Nx,3);
thresh = [hR + 0.02, phiR*XR + 0.02, phiR*(1-XR) + 0.02]; %threshold for shock detection

Y = Y0; %remove this later

%density_factor = phiL*XL*A.rhos1 + phiL*(1-XL)*A.rhos2 + 1;
%nul_eff_1 = A.nul*(1-phiL*XL/A.phimax)^(-2)/density_factor;

% we solve 0 = h_t + C_wm*(h^3/3)_x and 0 = n_t + (C_wm h^3/3 n)_x,
% where the system is non-dim. by the scales for the full system.
% Here C_wm is the ratio of kinematic viscosities of the fluid to the
% fluid mixture.
C_wm = (phiL*XL*A.rhos1 + phiL*(1-XL)*A.rhos2 + 1)/(1-phiL/A.phimax)^(-2);


%Find the appropriate (monodisperse) fluxes

%effective kinematic viscosity
nul_eff_1 = A.nul*(1-phiL*XL/A.phimax)^(-2)/(phiL*XL*A.rhos1 + 1); 

%effective fluid density (mixture of lighter particles and fluid)
rho_fluid_eff_1 = (1-phiL*XL)*A.rholiq + phiL*XL*A.rhopart1;

%Compute the relevant ODE constants
if(x_trans_2 > x_trans_1)
    A2 = set_constants2(ang,A.rhopart1,A.rhopart2,rho_fluid_eff_1);
    A2.nul = nul_eff_1; %not actually used; added for completeness
    A2.phi_shift = phiL*XL; %adjustment for viscosity functions, etc. in ODE
    
    %TODO: This setup doesn't make sense because the packing fraction, etc.
    %change so that one cannot just substitute a 'mixed' fluid into the sys.
    
    A2.ODEsolve = @ode45;
    [phi0s, fluxes_1, ~] = build_fluxtable1(A2,0,200);
    mono_flux = @(Y) bF1_s(Y,phi0s,fluxes_1); %monodisperse flux function
end

t = 0;
it = 1;

dt_max = tf/100;
set(0,'currentFigure',1);

%t_dim_1 = t_dim*A.nul/nul_eff_1;

C_wm_1 = (phiL*XL*A.rhos1 + 1)/(1-XL*phiL/A.phimax)^(-2);
%C_wm_1 = (phiL*XL*A.rhos1 + 1);
phase = 1;
max_speed = 1;
t_trans_1 = 0; t_trans_2 = 0;
% 
% h1l = [0];
% h1r = [0];
iout = 1;
nout = length(t_out);
Yout = zeros(length(X),length(Y(1,:)),length(t_out));
output_frame = false;
if(t_out<=t)
    Yout(:,:,iout) = Y(:,:);
    iout = iout + 1;
end

while(t < tf - eps)
    
    if(phase==1 && (fronts(it,1) >= x_trans_1))
       fprintf('\n **********End of transient 1************\n');
       t_trans_1 = t*t_dim;
       phase = 2;
    end
    
    if(phase==2 && (fronts(it,1) >= x_trans_2))
        t_trans_2 = t*t_dim;
        phase = 3;
       fprintf('\n ***********End of transient 2************\n');
    end
    
    max_speed = 0;
    if(phase==1) %Well mixed (fluid + p1 + p2)
        for i=1:Nx
            fluxes(i,1) = C_wm*Y(i,1)^3/3; %WELL MIXED FLUXES
            fluxes(i,2) = C_wm*Y(i,2)*Y(i,1)^2/3;
            fluxes(i,3) = C_wm*Y(i,3)*Y(i,1)^2/3;
            max_speed = max([max_speed,3*fluxes(i,:)/Y(i,1)]); %not quite a good estimate
        end
    elseif(phase==2) %p2 settled (fluid + p1 mixed)
        for i=1:Nx
            mflux = mono_flux(Y(i,[1 3]));
            fluxes(i,[1 3]) = C_wm_1*Y(i,1)^3*mflux; %fluid + heavier particles!
            fluxes(i,2) = C_wm_1*Y(i,2)*Y(i,1)^2*mflux(1); %lighter particles
            max_speed = max([max_speed,3*fluxes(i,:)/Y(i,1)]);
        end
    else %p1 and p2 settled
        for i=1:Nx
            fluxes(i,:) = bF(Y(i,:)',fpars{:});
            max_speed = max([max_speed,3*fluxes(i,:)/Y(i,1)]);
        end
    end
    
    %max_speed = 3*max(max(abs(fluxes./Y))); %estimate by max of f/h, g1/n1, g2/n2
    
    dt = CFL_max*dx/max_speed; %new dt
    dt = min(dt,dt_max);
    
    if(iout <= nout && (t > t_out(iout) - dt)) %Override if we want to integrate to an output frame
        dt = t_out(iout) - t;
        output_frame = 1;
    end
    
    a = dt/dx;
    
    set(0,'currentFigure',1);
    plot(X,Y);
    title(sprintf('t = %5.2f s',t*t_dim));
    drawnow
    
    %Y(1,:) = Y(1,:); %no change on left endpoint
    Y(2:end,:) = Y(2:end,:) - a*(fluxes(2:end,:) - fluxes(1:(end-1),:));
    %Y(2:(end-1),:) = (Y(1:(end-2),:) + Y(3:end,:))/2 ...
    %    - 0.5*a*(fluxes(3:end,:) - fluxes(1:(end-2),:));
    
    
    t = t + dt;
    
    if(output_frame)
        Yout(:,:,iout) = Y(:,:);
        fprintf('Output at t=%f, iout = %f\n',t,iout);
        iout = iout + 1;
        output_frame = false;
    end
 
    it = it + 1; 
    T(it) = t;
    if(mod(it,20)==1)
        fprintf('%f and dt = %f\n',t,dt);
    end
    
    %Compute the new front position
    for i=1:3
        itmp = find(Y(:,i) > thresh(i),1,'last');
        if(~isempty(itmp))
            fronts(it,i) = X(find(Y(:,i) > thresh(i),1,'last'));
        else
            fronts(it,i) = 0;
        end
    end
    
    %dt = 1e-4;
end

fprintf('Transient times: First: %f, Second: %f\n',t_trans_1,t_trans_2);


%Copy the rest of the results into the solution struct
bsol.T = T;
bsol.ang = ang;
bsol.t_out = t_out;
bsol.Y = Yout;
bsol.fronts = fronts;
%bsol.ictype = ic_type;
bsol.idata = [phiL XL hL; phiR XR hR];
bsol.Vsl = Vsl;
bsol.wtr = wtr;
bsol.res_length = res_length;
bsol.L = Ltr;
bsol.t_dim = t_dim;
% % bsol.h1r = h1r;
% % bsol.h1l = h1l;
%bsol.L_settled = L_settle;
%bsol.xf = xf;


% %Simple plot to show results (use a different script for better plots)
% figure
% set(0,'defaultTextInterpreter','tex')
% %Ytf = data_to_hp(Y(:,:,end));
% T = T*t_dim;
% plot(T,fronts(:,1),'-k',T,fronts(:,2),'-r',T,fronts(:,3),'-b');
% xlabel('T');
% legend('h front','n1','n2');
% %title(sprintf('Profile at t = %g, phiL = %g, XL = %g',tf,phiL,XL));


%------------Save the data----------------

if(~isempty(run_name))
    save(['./runs_autosave/',run_name,'.mat'],'bsol');
end