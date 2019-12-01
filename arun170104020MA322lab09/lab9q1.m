clear all;
clc;

a = 0;
b = 1;

%Part a
h = 0.01;
n = (b-a)/h;
x = a+h:h:b;
initial_value = 1;
final_value = 3;
p = @(x) (-x-1);
q = @(x) cos(x);
r = @(x) -exp(x);
central_sol_a = central_BVP(a, h, n, p, q, r, x, initial_value, final_value, 1);
forward_sol_a = forward_BVP(a, h, n, p, q, r, x, initial_value, final_value, 2);
backward_sol_a = backward_BVP(a, h, n, p, q, r, x, initial_value, final_value, 3);

x = [a x];
fprintf('Part a:\n');

fprintf('By central difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(central_sol_a)
    if mod(i,5) == 0
        fprintf('%d \t %f \t %f\n', i, x(i), central_sol_a(i));
    end
end

fprintf('By forward difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(forward_sol_a)
    if mod(i,5) == 0
        fprintf('%d \t %f \t %f\n', i, x(i), forward_sol_a(i));
    end
end

fprintf('By backward difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(backward_sol_a)
    if mod(i,5) == 0
     fprintf('%d \t %f \t %f\n', i, x(i), backward_sol_a(i));
    end
end

%Part b
h = 0.01;
n = (b-a)/h;
x = a+h:h:b;
initial_value = 0;
final_value = 0;
p = @(x) -2;
q = @(x) -1;
r = @(x) -x;
central_sol_b = central_BVP(a, h, n, p, q, r, x, initial_value, final_value, 4);
forward_sol_b = forward_BVP(a, h, n, p, q, r, x, initial_value, final_value, 5);
backward_sol_b = backward_BVP(a, h, n, p, q, r, x, initial_value, final_value, 6);

x = [a x];
fprintf('Part b:\n');

fprintf('By central difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(central_sol_b)
    if mod(i,5) == 0
     fprintf('%d \t %f \t %f\n', i, x(i), central_sol_b(i));
    end
end

fprintf('By forward difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(forward_sol_b)
    if mod(i,5) == 0
        fprintf('%d \t %f \t %f\n', i, x(i), forward_sol_b(i));
    end
end

fprintf('By backward difference:\n');
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(backward_sol_b)
    if mod(i,5) == 0
    fprintf('%d \t %f \t %f\n', i, x(i), backward_sol_b(i));
    end
end

function y = central_BVP(a, h, n, p, q, r, x, initial_value, final_value, fig)
    b(1) = r(x(1)) + ((1/h^2) + p(x(1))/(2*h))*initial_value;
    for i = 2:n-2
        b(i) = r(x(i));
    end
    b(n-1) = r(x(n-1)) - ((-1/h^2) + p(x(n-1))/(2*h))*final_value;
    %A = zeros(n-1,n-1);
    A(1,1) = (2/h^2) + q(x(1));
    A(1,2) = (-1/h^2) + p(x(1))/(2*h);
    for i = 2:n-2
        A(i,i-1) = (-1/h^2) - p(x(i))/(2*h);
        A(i,i) = (2/h^2) + q(x(i));
        A(i,i+1) = (-1/h^2) + p(x(i))/(2*h);
    end
    A(n-1,n-2) = (-1/h^2) - p(x(n-1))/(2*h);
    A(n-1,n-1) = (2/h^2) + q(x(n-1));
    b = b';
    y = A\b;
    y = y';
    y = [initial_value y final_value]; 
    y = y';
    figure(fig);
    x = [a x];
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
    title('Central difference');
end

function y = forward_BVP(a, h, n, p, q, r, x, initial_value, final_value, fig)
    b(1) = r(x(1)) + (1/h^2)*initial_value;
    for i = 2:n-2
        b(i) = r(x(i));
    end
    b(n-1) = r(x(n-1)) - ((-1/h^2) + p(x(n-1))/(h))*final_value;
    %A = zeros(n-1,n-1);
    A(1,1) = (2/h^2) - (p(x(1))/h) + q(x(1));
    A(1,2) = (-1/h^2) + p(x(1))/(h);
    for i = 2:n-2
        A(i,i-1) = (-1/h^2);
        A(i,i) = (2/h^2) + q(x(i)) - p(x(i))/h;
        A(i,i+1) = (-1/h^2) + p(x(i))/(h);
    end
    A(n-1,n-2) = (-1/h^2);
    A(n-1,n-1) = (2/h^2) + q(x(n-1)) - p(x(n-1))/h;
    b = b';
    y = A\b;
    y = y';
    y = [initial_value y final_value]; 
    y = y';
    figure(fig);
    x = [a x];
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
    title('Forward difference');
end

function y = backward_BVP(a, h, n, p, q, r, x, initial_value, final_value, fig)
    b(1) = r(x(1)) - (-(1/h^2) - p(x(1))/(h))*initial_value;
    for i = 2:n-2
        b(i) = r(x(i));
    end
    b(n-1) = r(x(n-1)) + (1/h^2)*final_value;
    %A = zeros(n-1,n-1);
    A(1,1) = (2/h^2) + q(x(1)) + p(x(1))/h;
    A(1,2) = (-1/h^2);
    for i = 2:n-2
        A(i,i-1) = (-1/h^2) - p(x(i))/(h);
        A(i,i) = (2/h^2) + q(x(i)) + p(x(i))/h;
        A(i,i+1) = (-1/h^2);
    end
    A(n-1,n-2) = (-1/h^2) - p(x(n-1))/(h);
    A(n-1,n-1) = (2/h^2) + q(x(n-1)) + p(x(n-1))/h;
    b = b';
    y = A\b;
    y = y';
    y = [initial_value y final_value]; 
    y = y';
    figure(fig);
    x = [a x];
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
    title('Backward difference');
end