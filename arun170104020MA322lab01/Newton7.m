function result = Newton(x)

while(abs(func(x))>0.000001)
    x = x - func(x)/derfunc(x);
end

disp('Solution by Newton method is  ' + string(x));

end

function f = func(x)
    f = 230* power(x,4) + 18*power(x,3) + 9*power(x,2) -221*x -9; 
end

function f = derfunc(x)
    f = 920*power(x,3) + 54*x*x + 18*x -221;
end


