%Simple bisection routine for the single-species shooting

function [x,it] = bisect(F,xL,xR,tol)

fxR = F(xR);
fxL = F(xL);

if(fxR*fxL > 0)
    fprintf('Warning: no bracketed root!\n');
    x = 0;
    it = 0;
    return
end

it = 0;
while(0.5*abs(xR - xL) > tol)
    it = it + 1;
    x = xL + 0.5*(xR - xL);
    fx = F(x);
    if(fxL*fx > 0)
        xL = x;
        fxL = fx;
    else
        xR = x;
    end
end