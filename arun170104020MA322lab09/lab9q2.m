clear all;
clc;

%Part a
a = -1;
b = 0;
h = 0.01;
n = (b-a)/h;
x = a:h:b;
initial_value = -2;
final_value = 1;
p = @(x) (-x);
q = @(x) 2;
r = @(x) -2 - (2 + x^2).*exp(x);
exact_solution = @(x) (x.^2 + x.*exp(x));
solution_a = Neumann_BVP(h, n, p, q, r, x, initial_value, final_value, exact_solution, 1);

%Part b
a = 0;
b = 1;
h = 0.01;
n = (b-a)/h;
x = a:h:b;
initial_value = -1;
final_value = 2*sin(1);
p = @(x) 0;
q = @(x) -x;
r = @(x) -(3 - x - x.^2 + x.^3).*sin(x) - (4*x).*cos(x);
exact_solution = @(x) (sin(pi*x)).^2;
solution_b = Without_exact_sol_Neumann_BVP(h, n, p, q, r, x, initial_value, final_value, 3);

function y = Neumann_BVP(h, n, p, q, r, x, initial_value, final_value, exact_solution, fig)  
    b(1) = initial_value*h;
    for i = 2:n
        b(i) = r(x(i));
    end
    b(n+1) = final_value*h;
    A(1,1) = -1;
    A(1,2) = 1;
    for i = 2:n
        A(i,i-1) = (-1/h^2) - p(x(i))/(2*h);
        A(i,i) = (2/h^2) + q(x(i));
        A(i,i+1) = (-1/h^2) + p(x(i))/(2*h);
    end
    A(n+1,n) = -1;
    A(n+1,n+1) = 1;
    b = b';
    y = A\b;
    err = abs(y - exact_solution(x)');
    figure(fig);
    plot(x,y);
    hold on;
    plot(x,exact_solution(x),'--', 'linewidth', 2);
    legend('Approximation', 'Exact Solution');
    xlabel('x(i)');
    ylabel('y(x(i))');
    figure(fig+1);
    plot(x,err);
    xlabel('x(i)');
    ylabel('Error');
    fprintf('x \t\t Approx \t Exact \t\t Error\n');
    for i = 1:length(y)
        if mod(i,5) == 0
         fprintf('%f \t %f \t %f \t %e\n', x(i), y(i), exact_solution(x(i)), err(i));
        end
    end
end

function y = Without_exact_sol_Neumann_BVP(h, n, p, q, r, x, initial_value, final_value, fig)  
    b(1) = initial_value*h;
    for i = 2:n
        b(i) = r(x(i));
    end
    b(n+1) = final_value*h;
    A(1,1) = -1;
    A(1,2) = 1;
    for i = 2:n
        A(i,i-1) = (-1/h^2) - p(x(i))/(2*h);
        A(i,i) = (2/h^2) + q(x(i));
        A(i,i+1) = (-1/h^2) + p(x(i))/(2*h);
    end
    A(n+1,n) = -1;
    A(n+1,n+1) = 1;
    b = b';
    y = A\b;
    figure(fig);
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
    fprintf('\nx \t\t Approx\n');
    for i = 1:length(y)
        if mod(i,5) == 0
            fprintf('%f \t %f \n', x(i), y(i));
        end
    end
end