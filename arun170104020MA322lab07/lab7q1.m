clear all;
clc;
clf;

f = @(lmd,y,t)  (lmd*y + (1-lmd)*cos(t) - (1+lmd)*sin(t));
sol = @(t) (sin(t) + cos(t));

%Part a with lambda = -5, h = 0.01, 0.001, 0.005, 0.0025
% explicit_euler(-5,0.01,f,sol,1,0,10);
% explicit_euler(-5,0.001,f,sol,1,0,10);
% explicit_euler(-5,0.005,f,sol,1,0,10);
% explicit_euler(-5,0.0025,f,sol,1,0,10);
% implicit_euler(-5,0.01,f,sol,1,0,10);
% implicit_euler(-5,0.001,f,sol,1,0,10);
% implicit_euler(-5,0.005,f,sol,1,0,10);
% implicit_euler(-5,0.0025,f,sol,1,0,10);

%Part b with lambda = 5, h = 0.01,0.001
% explicit_euler(5,0.01,f,sol,1,0,10);
% explicit_euler(5,0.001,f,sol,1,0,10);
% implicit_euler(5,0.01,f,sol,1,0,10);
% implicit_euler(5,0.001,f,sol,1,0,10);

function explicit_euler(lmd,h,f,sol,init,a,b)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    for i = 1:n
        y(i) = prev + h*f(lmd,prev,t);
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
    figure(3);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end

function implicit_euler(lmd,h,f,sol,init,a,b)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    tol = 1e-8;
    for i=1:n
        y_old = prev + h*f(lmd,prev,t+h);
        y_new = prev + h*f(lmd,y_old,t+h);
        while abs(y_new - y_old) > tol
            save = y_new;
            y_new = prev + h*f(lmd,y_old,t+h);
            y_old = save;
        end
        y(i) = y_new;
        prev = y(i);
        t = t+h;
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
    figure(3);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end