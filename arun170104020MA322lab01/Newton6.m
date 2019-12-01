function result = Newton(x,y)
switch nargin
    case 1
        y = -1.5;
end

while(abs(func(x,y))>1e-8)
    x = x - func(x,y)/derfunc(x);
end

disp('Solution by Newton method is ' + string(x));

end

function f = func(x,y)
    f = exp(x) + y - atan(x); 
end

function f = derfunc(x)
    f = exp(x) - (1/(1+x*x));
end
