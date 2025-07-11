%%Sets constants for the bidisperse system.
% Input angle is in degrees.

%TODO: make varargin take kwargs!

function A = set_constants2(varargin)

%ANGLE rhop1 rhop2 rholiq
ang = varargin{1};
A.ang = ang;
phimax = 0.61;  %Maximum packing fraction

alpha = pi * ang/180; % inclination angle
Kcoll = 0.41;
%Common diameters
%GSB3: 0.725mm
%GSB5: 0.3625mm
%GSB7: 0.215mm 

if(nargin > 1 && ~isempty(varargin{2}))
    rhopart1 = varargin{2};  % density of smaller particle
    rhopart2 = varargin{3}; % density of larger particle
    rholiq = varargin{4}; 
else   
    d1 = 625 * 10^-6;  % diameter of larger particle (m) (new particles)
    d2 = 200 * 10^-6; % diameter of smaller particle (GSB5)
    rholiq = 971;    % density of the liquid (kg/m^3)
    rho = 2475; %density of the sphere
    rhopart1 = rho;
    rhopart2 = rho;
end

if(nargin > 4)
    A.phi_shift = varargin{5};
else
    A.phi_shift = 0;
end


    

rhos = (rho-rholiq)/rholiq; %Normalized density
    
Kvisc=0.62; % shear induced migration constants

%Tracer diffusivity
beta = 2; Kt = 0.5;
A.nul = 1e-3*971/rholiq; % kinematic viscosity in m^2/s (1000 cSt)

%Dilute coefficient (Acrivos 1987)
%A.Dtr = @(phi) Kt*phi.^beta;

%Stokesian prediction (Sierou and Brady 2004)
phia = 0.4; %cutoff after which Dtr is const.
A.Dtr = @(phi) (phi <= phia)*Kt*phi.^2 + (phi > phia)*Kt*(phia)^2;


c1   = 2*(Kvisc-Kcoll)/Kcoll;      % constant c1 in ode
c2   = 2/(9*Kcoll)*cot(alpha);     % constant c2 in ode
%A.c0 = 2*(rhos2-rhos1)*cot(alpha)/9; %constant c0 in X' ODE (two species)
A.c0 = 2*cot(alpha)/9; %%%%%%%THIS IS A PLACEHOLDER VALUE%%%%%%%%%%


%Critical phi as a function of density
A.phic = @(rhos) -(c2*rhos+1)./(2*rhos)+sqrt(c2+(0.5*(c2*rhos+1)./rhos).^2);

%Order of pole for 1/(phimax - phi) at z=1 for ridged case
A.p_r = @(rhos) (phimax + rhos*phimax^2 + c2*rhos*(phimax-1))/(c1*phimax*(1+rhos*phimax));

%IVP solver settings
A.vopt = odeset('RelTol', 1e-7, 'AbsTol', 5e-8,'MaxStep', 5e-2);
A.ODEsolve = @ode15s;
A.gamma = A.c0/Kt; %Scaling factor
A.max_it = 100;



A.rhopart1 = rhopart1; A.rhopart2 = rhopart2; A.rholiq = rholiq;
A.phimax = phimax; 
A.alpha = alpha;
A.d1 = d1; A.d2 = d2;
A.c1 = c1; A.c2 = c2;
A.rho = rho; 
A.rhos = rhos;
A.rhos1 = rhos;
A.rhos2 = rhos; %These are all the same for bidisperse case (all the same
%density)
A.Kv = Kvisc; A.Kc = Kcoll; A.Kt = Kt; A.beta = beta;