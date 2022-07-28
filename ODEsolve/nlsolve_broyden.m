%Broyden's method with a partially implemented line search (quadratic only)
%The solver works with a good initial guess, but will sometimes
%fail if the guess is not close or it otherwise gets confused.

%Eventually, this should be replaced by a better non-linear solver,
%since it can be a bit unreliable.

% Solves the system F(x) = 0 using Broyden's method with a partially
% implemented quadratic linesearch. Enforces bounds on the domain in a 
% simple (but crude) way. 
%
% TODO: Update QR instead of the actual Jacobian [see Numerical Recipes]
%
% Input:
%   x0       -  The initial guess (column vector)
%   F        -  The given total/lighter concentration
%   n        -  The dimension of the system
%   tol      -  Stopping tolerance
%   bounds   -  nx2 array specifying the bounds of the domain.
%               This isn't handled gracefully, but can sometimes stop
%               the solver from choosing 'bad' points.
%   max_it   -  max. number of iterations
%   pars -  (optional) parameters for F, as a comma-separated list (if
%               none, set pars = {} )
%
% Output:
%   x        -  the computed solution
%   it       -  the number of iterations

function [x,it]=nlsolve_broyden(x0,F,n,tol,bounds,max_it,pars)

verb = 0;   %1 = show final point and iterations
            %2 = show point at each itereation
EPS = 1e-7;
cond_min = 1e-9; %minimum inverse condition number (solver resets in this case)
tf = 0.5; %Used for the linesearch (backtracking factor)
Jac = eye(n); %Initial Jacobian approximation

x = x0; Fx = F(x,pars{:});
it = 0; 

needs_reset = 1;

%Make sure that we only search within specified bounds.
if(~isempty(bounds))
    out_of_bounds = @(x) any(x - bounds(:,1) <= EPS) | any(bounds(:,2) - x <= EPS);
else
    out_of_bounds = @(x) 0;
end
    
if(out_of_bounds(x0))
    fprintf('x0 is not within bounds...\n');
    x = (bounds(:,1) + bounds(:,2))/2;
end

%Main loop
while(norm(Fx)>tol && it < max_it)
    it = it + 1;
    if(verb>=2)
        fprintf('it = %i\t xn = %f\t%f\t',it,x(1),x(2));
        fprintf('g1: %f \t g2: %f\n',Fx(1),Fx(2));
    end
    if(mod(it,30)==0) %arbitrary time to give up
        %If no good progress, recompute Jacobian
        needs_reset = 1;
    end
    
    if(rcond(Jac) < cond_min || needs_reset)
        needs_reset = 0;
        Jac = eye(n);
    end
    
    del = Jac\(-Fx); %Since n=2, this is not too bad
    tau = 1;
    xn = x + tau*del;
    
    %Cheap way to enforce bounds.
    while( out_of_bounds(xn) )
        tau = tf*tau;
        xn = x + tau*del;
    end    
    del = tau*del;
    
    Fx_old = Fx;
    Fx = F(xn,pars{:});
    if(isnan(norm(Fx)))
        fprintf('Solver has failed.  (Broyden)\n');
        it = -1;
        x = zeros(n,1);
        return
    end
    if(norm(Fx) > norm(Fx_old))
        %Do a quadratic linesearch here to get a better value
        c0 = norm(Fx_old)^2; %g(0)
        c1 = 2*dot(Fx_old,Jac*del); %g'(0)
        c2 = norm(Fx)^2 - c1 - c0;      
        if(c2~=0 && c1~=0)
            tau = -c1/(2*c2);
        else
           tau = 0.5;
        end
        tau = min(tau,1);
        tau = max(tau,0.1); %Don't take too small a step
    end
    del = tau*del; 
    
    x = x + del; %update x
    Fx = F(x,pars{:}); 
    v = Fx - Fx_old - Jac*del;
    nd = norm(del)^2;
    if(nd==0)
        needs_reset = 1;
    else
        for i=1:n
            for j=1:n
                Jac(i,j) = Jac(i,j) + v(i)*del(j)/nd;
            end
        end
    end
end

if(verb>=1)
    fprintf('x0 = %f %f, it = %i\n',x0(1),x0(2),it);
end
    