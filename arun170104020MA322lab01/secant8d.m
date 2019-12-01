function result = secant(x,y)

i=2;
xn(1) = x;
xn(2) = y;
while(abs(func(xn(i)))>1e-4)
    i = i+1;
    xn(i) = xn(i-1) - (func(xn(i-1))*(xn(i-1)-xn(i-2)))/(func(xn(i-1))-func(xn(i-2)));
end

disp('Solution by Secant method is  ' + string(xn(i)));
end

function f = func(x)
    f = x*x*x + 4.001*x*x + 4.002*x + 1.101; 
end
