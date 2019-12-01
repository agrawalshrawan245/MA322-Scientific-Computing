clear all;
clc;
clf;

a = 0;
b = 1;

%Part a
h = 0.01;
n = (b-a)/h;
x = a:h:b;
initial_value = 1;
final_value = 1;
p = @(x) 0;
q = @(x) x;
r = @(x) -1;
solution_1 = mixed_BVP_1(h, n, p, q, r, x, initial_value, final_value, 1);
fprintf('n \t x(i) \t\t Approximation\n');
for i =1:length(solution_1)
    if mod(i,5) == 0
     fprintf('%d \t %f \t %f\n', i, x(i), solution_1(i));
    end
end

%Part b
h = 0.01;
n = (b-a)/h;
x = a:h:b;
initial_value = 1;
final_value = 1;
p = @(x) 3;
q = @(x) -2;
r = @(x) -2;
solution_2 = mixed_BVP_2(h, n, p, q, r, x, initial_value, final_value, 2);
fprintf('\nn \t x(i) \t\t Approximation\n');
for i =1:length(solution_2)
    if mod(i,5) == 0
        fprintf('%d \t %f \t %f\n', i, x(i), solution_2(i));
    end
end

function y = mixed_BVP_1(h, n, p, q, r, x, initial_value, final_value, fig)  
    b(1) = initial_value;
    for i = 2:n
        b(i) = r(x(i));
    end
    b(n+1) = final_value;
    A(1,1) = 1-(1/h);
    A(1,2) = 1/h;
    for i = 2:n
        A(i,i-1) = (-1/h^2) - p(x(i))/(2*h);
        A(i,i) = (2/h^2) + q(x(i));
        A(i,i+1) = (-1/h^2) + p(x(i))/(2*h);
    end
    A(n+1,n+1) = 1;
    b = b';
    y = Gauss_Siedel(A,b);
    figure(fig);
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
end

function y = mixed_BVP_2(h, n, p, q, r, x, initial_value, final_value, fig)  
    b(1) = initial_value;
    for i = 2:n
        b(i) = r(x(i));
    end
    b(n+1) = final_value;
    A(1,1) = -1-(1/h);
    A(1,2) = 1/h;
    for i = 2:n
        A(i,i-1) = (-1/h^2) - p(x(i))/(2*h);
        A(i,i) = (2/h^2) + q(x(i));
        A(i,i+1) = (-1/h^2) + p(x(i))/(2*h);
    end
    A(n+1,n) = (-1/h);
    A(n+1,n+1) = 1 + (1/h);
    b = b';
    y = A\b;
    figure(fig);
    plot(x,y);
    xlabel('x(i)');
    ylabel('y(x(i))');
end

function x = Gauss_Siedel(A, b)
    x=zeros(length(b),1);
    n=size(x,1);
    normVal=Inf; 
    tol=1e-8;
    while normVal>tol
        x_old=x;
        for i=1:n
            sigma=0;
            for j=1:i-1
                    sigma=sigma+A(i,j)*x(j);
            end
            for j=i+1:n
                    sigma=sigma+A(i,j)*x_old(j);
            end
            x(i)=(1/A(i,i))*(b(i)-sigma);
        end
        normVal=norm(x_old-x);
    end
end