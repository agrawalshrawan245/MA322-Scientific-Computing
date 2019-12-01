clear all;
clc;
clf;

%Given P.D.E. as ut = alpha*uxx

%Part a
exact_sol = @(x,t) (exp(-t).*sin((pi/2)*x) + exp(-t/4).*sin((pi/4)*x));
initial_cond = @(x) (sin((pi/4)*x).*(1 + 2*cos((pi/4)*x)));
h = 0.2;
k = 0.02;
xl = 0;
xr = 4;
x = xl:h:xr;
ti = 0;
tf = 1;
t = ti:k:tf;
bcl = @(t) 0;
bcr = @(t) 0;
alpha = (4/(pi*pi));
% FTCS_final_time_sol_a = FTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% BTCS_final_time_sol_a = BTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Crank_Nicolson_final_time_sol_a = Crank_Nicolson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Richardson_final_time_sol_a = Richardson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Du_Fort_Frankel_final_time_sol_a = Du_Fort_Frankel(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);

%Part b
exact_sol = @(x,t) (exp(-(pi*pi)*t).*sin(pi*x));
initial_cond = @(x) (sin(pi*x));
h = 0.1;
k = 0.002;
xl = 0;
xr = 1;
x = xl:h:xr;
ti = 0;
tf = 1;
t = ti:k:tf;
bcl = @(t) 0;
bcr = @(t) 0;
alpha = 1;
% FTCS_final_time_sol_b = FTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% BTCS_final_time_sol_b = BTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Crank_Nicolson_final_time_sol_b = Crank_Nicolson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Richardson_final_time_sol_b = Richardson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Du_Fort_Frankel_final_time_sol_b = Du_Fort_Frankel(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);

%Part c
exact_sol = @(x,t) (exp(-((pi*pi*t)/4)).*sin((pi*x)/2) + (0.5)*(exp(-4*pi*pi*t).*sin(2*pi*x)));
initial_cond = @(x) (sin((pi*x)/2) + (0.5)*sin(2*pi*x));
h = 0.1;
k = 0.002;
xl = 0;
xr = 1;
x = xl:h:xr;
ti = 0;
tf = 1;
t = ti:k:tf;
bcl = @(t) 0;
bcr = @(t) exp(-((pi*pi*t)/4));
alpha = 1;
% FTCS_final_time_sol_c = FTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% BTCS_final_time_sol_c = BTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Crank_Nicolson_final_time_sol_c = Crank_Nicolson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
 Richardson_final_time_sol_c = Richardson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);
% Du_Fort_Frankel_final_time_sol_c = Du_Fort_Frankel(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, 1);

function final_time_sol = FTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, fig)
    n = (tf-ti)/k;
    max_err = [];
    for m = 2:10
        if m>2 
            clear u;
        end
        h = (xr - xl)/m;
        lambda = (k/h^2)*alpha;
        u(1,:) = initial_cond(xl:h:xr); 
        for i = 2:n+1
                u(i,1) = bcl(ti + (i-1)*k);
                u(i, m+1) = bcr(ti + (i-1)*k);
                for j = 2:m
                        u(i,j) = lambda*u(i-1,j-1) + (1 - 2*lambda)*u(i-1,j) + lambda*u(i-1,j+1);
                end
        end
        [X,T] = meshgrid(xl:h:xr,ti:k:tf);
        err = max(max(abs(u - exact_sol(X,T))));
        max_err = [max_err err];
    end
    x = xl:h:xr;
    final_time_sol = u(n+1,:);
    figure(fig);
    plot(x,final_time_sol);
    hold on;
    plot(x,exact_sol(x,tf),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('u(x(i))');
    figure(fig+1);
    [X,T] = meshgrid(x,t);
    Z = exact_sol(X,T);
    surf(X,T,Z);
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    title('Surface plot of exact solution');
    figure(fig+2);
    [X,T] = meshgrid(x,t);
    surf(X,T,u);
    title('Surface plot of numerical solution');
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    figure(fig+3);
    loglog(2:10,max_err);
    xlabel('N');
    ylabel('Max error');
end

function final_time_sol = BTCS(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, fig)
   n = (tf-ti)/k;
   max_err = [];
    for m = 3:20
        if m>3 
            clear u;
            clear solution_matrix;
        end
        h = (xr - xl)/m;
        lambda = (k/h^2)*alpha;
        A0(1,1) = 1 + 2*lambda;
        A0(1,2) = -lambda;
        for i = 2:m-2
            A0(i,i-1) = -lambda;
            A0(i,i) = 1 + 2*lambda;
            A0(i,i+1) = -lambda;
        end
        A0(m-1,m-2) = -lambda;
        A0(m-1,m-1) = 1 + 2*lambda;
        u(:,1) = initial_cond(xl+h:h:xr-h);
        for i = 1:n
            u(1,i) = u(1,i) + lambda*bcl(ti + i*k);
            u(m-1,i) = u(m-1,i) + lambda*bcr(ti + i*k);
            u(:,i+1) = A0\u(:,i);
            l = bcl(ti + (i-1)*k);
            r = bcr(ti + (i-1)*k);
            solution_matrix(:,i) = vertcat(l, u(:,i), r);
        end
        l = bcl(tf);
        r = bcr(tf);
        solution_matrix(:,n+1) = vertcat(l, u(:,n+1), r);
        [X,T] = meshgrid(xl:h:xr,ti:k:tf);
        err = max(max(abs(solution_matrix' - exact_sol(X,T))));
        max_err = [max_err err];
    end
    x = xl:h:xr;
    final_time_sol = solution_matrix(:,n+1);
    figure(fig);
    plot(x,final_time_sol);
    hold on;
    plot(x,exact_sol(x,tf),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('u(x(i))');
    figure(fig+1);
    [X,T] = meshgrid(x,t);
    Z = exact_sol(X,T);
    surf(X,T,Z);
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    title('Surface plot of exact solution');
    figure(fig+2);
    [X,T] = meshgrid(x,t);
    Z = solution_matrix';
    surf(X,T,Z);
    title('Surface plot of numerical solution');
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    figure(fig+3);
    loglog(3:20,max_err);
    xlabel('N');
    ylabel('Max error');
end

function final_time_sol = Crank_Nicolson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, fig)
   n = (tf-ti)/k;
   max_err = [];
    for m = 3:20
        if m>3 
            clear u;
            clear solution_matrix;
        end
        h = (xr - xl)/m;
            lambda = (k/h^2)*alpha;
        A0(1,1) = 1 + lambda;
        A0(1,2) = -lambda/2;
        for i = 2:m-2
            A0(i,i-1) = -lambda/2;
            A0(i,i) = 1 + lambda;
            A0(i,i+1) = -lambda/2;
        end
        A0(m-1,m-2) = -lambda/2;
        A0(m-1,m-1) = 1 + lambda;
        A1(1,1) = 1 - lambda;
        A1(1,2) = lambda/2;
        for i = 2:m-2
            A1(i,i-1) = lambda/2;
            A1(i,i) = 1 - lambda;
            A1(i,i+1) = lambda/2;
        end
        A1(m-1,m-2) = lambda/2;
        A1(m-1,m-1) = 1 - lambda;
        u(:,1) = initial_cond(xl+h:h:xr-h);
        for i = 1:n
            b = A1*u(:,i);
            b(1) = b(1) + (lambda/2)*bcl(ti + (i-1)*k) + (lambda/2)*bcl(ti + i*k);
            b(m-1) = b(m-1) + (lambda/2)*bcr(ti + (i-1)*k) + (lambda/2)*bcr(ti + i*k);
            u(:,i+1) = A0\b;
            l = bcl(ti + (i-1)*k);
            r = bcr(ti + (i-1)*k);
            solution_matrix(:,i) = vertcat(l, u(:,i), r);
        end
        l = bcl(tf);
        r = bcr(tf);
        solution_matrix(:,n+1) = vertcat(l, u(:,n+1), r);
        [X,T] = meshgrid(xl:h:xr,ti:k:tf);
        err = max(max(abs(solution_matrix' - exact_sol(X,T))));
        max_err = [max_err err];
    end
    x = xl:h:xr;
    final_time_sol = solution_matrix(:,n+1);
    figure(fig);
    plot(x,final_time_sol);
    hold on;
    plot(x,exact_sol(x,tf),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('u(x(i))');
    figure(fig+1);
    [X,T] = meshgrid(x,t);
    Z = exact_sol(X,T);
    surf(X,T,Z);
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    title('Surface plot of exact solution');
    figure(fig+2);
    [X,T] = meshgrid(x,t);
    Z = solution_matrix';
    surf(X,T,Z);
    title('Surface plot of numerical solution');
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    figure(fig+3);
    loglog(3:20,max_err);
    xlabel('N');
    ylabel('Max error');
end

function final_time_sol = Richardson(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, fig)
    n = (tf-ti)/k;
    max_err = [];
    for m = 2:10
        if m>2 
            clear u;
        end
        h = (xr - xl)/m;
        lambda = (k/h^2)*alpha;
        u(1,:) = initial_cond(xl:h:xr);  
        for i = 2:n+1
            u(i,1) = bcl(ti + (i-1)*k);
            u(i, m+1) = bcr(ti + (i-1)*k);
            for j = 2:m
                if i==2
                    u(i,j) = lambda*u(i-1,j-1) + (1 - 2*lambda)*u(i-1,j) + lambda*u(i-1,j+1);
                else
                    u(i,j) = u(i-2,j) + 2*lambda*(u(i-1,j+1) - 2*u(i-1,j) + u(i-1, j-1));
                end
            end
        end
        [X,T] = meshgrid(xl:h:xr,ti:k:tf);
        err = max(max(abs(u - exact_sol(X,T))));
        max_err = [max_err err];
    end
    x = xl:h:xr;
    final_time_sol = u(n+1,:);
    figure(fig);
    plot(x,final_time_sol);
    hold on;
    plot(x,exact_sol(x,tf),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('u(x(i))');
    figure(fig+1);
    [X,T] = meshgrid(x,t);
    Z = exact_sol(X,T);
    surf(X,T,Z);
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    title('Surface plot of exact solution');
    figure(fig+2);
    [X,T] = meshgrid(x,t);
    Z = u;
    surf(X,T,Z);
    title('Surface plot of numerical solution');
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    figure(fig+3);
    plot(2:10,max_err);
    xlabel('N');
    ylabel('Max error');
end

function final_time_sol = Du_Fort_Frankel(x, t, h ,k, xl, xr, ti, tf, alpha, initial_cond, bcl, bcr, exact_sol, fig)
    n = (tf-ti)/k;
    max_err = [];
    for m = 2:20
        if m>2 
            clear u;
        end
        h = (xr - xl)/m;
        lambda = (k/h^2)*alpha;
        u(1,:) = initial_cond(xl:h:xr); 
        for i = 2:n+1
            u(i,1) = bcl(ti + (i-1)*k);
            u(i, m+1) = bcr(ti + (i-1)*k);
            for j = 2:m
                if i==2
                    u(i,j) = lambda*u(i-1,j-1) + (1 - 2*lambda)*u(i-1,j) + lambda*u(i-1,j+1);
                else
                    u(i,j) = u(i-2,j)*((1 - 2*lambda)/(1 + 2*lambda)) + ((2*lambda)/(1 + 2*lambda))*(u(i-1,j+1) + u(i-1,j-1));
                end
            end
        end
        [X,T] = meshgrid(xl:h:xr,ti:k:tf);
        err = max(max(abs(u - exact_sol(X,T))));
        max_err = [max_err err];
    end
    x = xl:h:xr;
    final_time_sol = u(n+1,:);
    figure(fig);
    plot(x,final_time_sol);
    hold on;
    plot(x,exact_sol(x,tf),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('u(x(i))');
    figure(fig+1);
    [X,T] = meshgrid(x,t);
    Z = exact_sol(X,T);
    surf(X,T,Z);
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    title('Surface plot of exact solution');
    figure(fig+2);
    [X,T] = meshgrid(x,t);
    Z = u;
    surf(X,T,Z);
    title('Surface plot of numerical solution');
    xlabel('x(i)');
    ylabel('t(i))');
    zlabel('u(x(i), t(i))');
    figure(fig+3);
    loglog(2:20,max_err);
    xlabel('N');
    ylabel('Max error');
end