clear all;
clc;
clf;

f = @(t,y) (-y + (1.1 + t)*t^(0.1));
sol = @(t) (t.^(1.1));
a = 0;
b = 5;

max_err(1) = RK2(0.05, f, sol, 0, a, b, 1);
max_err(2) = RK2(0.025, f, sol, 0, a, b, 3);
max_err(3) = RK2(0.0125, f, sol, 0, a, b, 5);
max_err(4) = RK2(0.00625, f, sol, 0, a, b, 7);

starting_h = 0.00625;
fprintf('h \t\t Max Error \t log2(E(n)/E(2n))\n');
for i = 1:4
    fprintf('%f \t %f', power(2,4-i)*starting_h, max_err(i));
    if(i<4)
        fprintf('\t %f',log2(max_err(i)/max_err(i+1)));
    end
    fprintf('\n');
end

function max_err = RK2(h,f,sol,init,a,b,fig)
    prev = init;
    n = (b-a)/h;
    fprintf('For h = %.5f:\n', h);
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    for i = 1:n
        y(i) = prev + (h/4)*(f(t, prev) + 3*f(t + (2/3)*h, prev + (2/3)*h*f(t, prev)));
        prev = y(i);
        t = t + h;
        err(i) = abs(sol(t) - y(i));
        if(mod(i,(1/h))==0)
             fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));
        end  
    end    
    y = [init y];
    err = [0 err];
    figure(fig);
    plot(a:h:b, y);
    hold on;
    xx = a:h:b;
    plot(xx,sol(xx),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('t(i)');
    ylabel('y(t(i))');
    figure(fig+1);
    plot(a:h:b,err);
    xlabel('t(i)');
    ylabel('Error');
    max_err = max(err);
end