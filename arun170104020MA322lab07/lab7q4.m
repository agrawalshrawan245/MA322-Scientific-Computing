clear all;
clc;
clf;

% %Part a
f = @(t,y) (2 - 2*t.*y)./(1 + t.*t);
sol = @(t) (2.*t + 1)./(t.*t + 1);
a = 0;
b = 1;
init = 1;
h = 0.01;
explicit_euler(h,f,sol,init,a,b,1);
implicit_euler(h,f,sol,init,a,b,4);
modified_euler(h,f,sol,init,a,b,7);
trapezoidal(h,f,sol,init,a,b,10);

%Part b
% f = @(t,y) (y^2 + y)/t;
% sol = @(t) (2.*t)./(1-2.*t);
% a = 1;
% b = 3;
% init = -2;
% h = 0.02;
% explicit_euler(h,f,sol,init,a,b,1);
% implicit_euler(h,f,sol,init,a,b,4);
% modified_euler(h,f,sol,init,a,b,7);
% trapezoidal(h,f,sol,init,a,b,10);

%Part c
% f = @(t,y) 1 + (y/t) + (y/t)^2;
% sol = @(t) t.*tan(log(t));
% a = 1;
% b = 3;
% init = 0;
% h = 0.002;
% explicit_euler(h,f,sol,init,a,b,1);
% implicit_euler(h,f,sol,init,a,b,4);
% modified_euler(h,f,sol,init,a,b,7);
% trapezoidal(h,f,sol,init,a,b,10);

%Part d
% f = @(t,y) exp(t-y);
% sol = @(t) log(exp(t) + exp(1) -1);
% a = 0;
% b = 1;
% init = 1;
% h = 0.005;
% explicit_euler(h,f,sol,init,a,b,1);
% implicit_euler(h,f,sol,init,a,b,4);
% modified_euler(h,f,sol,init,a,b,7);
% trapezoidal(h,f,sol,init,a,b,10);

function explicit_euler(h,f,sol,init,a,b,fig)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    for i = 1:n
        y(i) = prev + h*f(t,prev);
        prev = y(i);
        t = t + h;
        err(i) = abs(sol(t) - y(i));
        fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));
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
    figure(fig+2);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end

function implicit_euler(h,f,sol,init,a,b,fig)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    tol = 1e-8;
    for i=1:n
        y_old = prev + h*f(t+h,prev);
        y_new = prev + h*f(t+h,y_old);
        while abs(y_new - y_old) > tol
            save = y_new;
            y_new = prev + h*f(t+h,y_old);
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
    figure(fig+2);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end

function modified_euler(h,f,sol,init,a,b,fig)
    prev = init;
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    for i = 1:n
        y(i) = prev + (h/2)*(f(t,prev) + f(t+h,prev+h*(f(t,prev))));
        prev = y(i);
        t = t + h;
        err(i) = abs(sol(t) - y(i));
        fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));
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
    figure(fig+2);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end

function trapezoidal(h,f,sol,init,a,b,fig)
    n = (b-a)/h;
    fprintf('t \t\t Approx \t Exact \t\t Error\n');
    t = a;
    prev = init;
    tol = 1e-8;
    for i=1:n
        y_old = prev + (h/2)*(f(t,prev)+f(t+h,prev));
        y_new = prev + (h/2)*(f(t,prev) + f(t+h,y_old));
        while abs(y_new - y_old) > tol
            save = y_new;
            y_new = prev + (h/2)*(f(t,prev) + f(t+h,y_old));
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
    figure(fig+2);
    semilogx(0:n,err);
    title('Semilogx plot of N vs. Error');
end