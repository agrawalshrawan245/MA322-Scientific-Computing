clear all;
clc;
clf;

f = @(t, y) (y/t + t*sec(y/t));
sol = @(t) (t.*asin(t));

a = 0;
b = 1;
h = power(2,-7);

RK4(h, f, sol, 0, a, b);

function RK4(h, f, sol, init, a, b)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    for i = 1:n
        if (t==0 && prev == 0)
            k1 = 0;
        else      
            k1 = h*f(t,prev);
        end
        k2 = h*f(t + (h/2), prev + (k1/2));
        k3 = h*f(t + (h/2), prev + (k2/2));
        k4 = h*f(t + h, prev + k3);
        y(i) = prev + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        prev = y(i);
        t = t + h;
        err(i) = abs(sol(t) - y(i));
        fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));  
    end    
    y = [init y];
    err = [0 err];
    figure(1);
    plot(a:h:b, y);
    hold on;
    xx = a:h:b;
    plot(xx,sol(xx),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('t(i)');
    ylabel('y(t(i))');
    figure(2);
    plot(a:h:b,err);
    xlabel('t(i)');
    ylabel('Error');
end