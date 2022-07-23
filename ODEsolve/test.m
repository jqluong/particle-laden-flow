f = @(x,y) x^2 + y^2;
g = @(x,y) f(x,y) * f(x,y);
g(1,1)