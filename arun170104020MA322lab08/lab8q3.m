clear all;
clc;

k = 6.22*power(10,-19);
n1 = 2*power(10,3);
n2 = 2*power(10,3);
n3 = 3*power(10,3);

f = @(t,y) k*((n1-(y/2)).^2)*((n2 - (y/2)).^2)*((n3 - ((3*y)/4)).^3);

a = 0;
b = 0.2;
h = 1e-3;

RK4(h, f, 0, a, b);

function RK4(h, f, init, a, b)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx\n');
    t = a;
    for i = 1:n    
        k1 = h*f(t,prev);
        k2 = h*f(t + (h/2), prev + (k1/2));
        k3 = h*f(t + (h/2), prev + (k2/2));
        k4 = h*f(t + h, prev + k3);
        y(i) = prev + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        prev = y(i);
        t = t + h;
        fprintf('%f \t %f \n', t,y(i));  
    end    
end