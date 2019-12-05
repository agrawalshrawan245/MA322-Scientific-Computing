clear all;
clc;

% f = @(x) [6*x(1) - 2*cos(x(2)*x(3)) - 1; 9*(x(2)) + sqrt(x(1)^2 + sin(x(3)) + 1.06) + 0.9; 60*x(3) + 3*exp(-x(1)*x(2)) + 10*pi - 3];
% Df = @(x) [6, 2*x(3)*sin(x(2)*x(3)), 2*x(2)*sin(x(2)*x(3));  x(1)/(sqrt(x(1)^2 + sin(x(3)) + 1.06)), 9, cos(x(3))/(sqrt(x(1)^2 + sin(x(3)) + 1.06)); -x(2)*3*exp(-x(1)*x(2)), -x(1)*3*exp(-x(1)*x(2)), 60];
% 
% x_old = [0 0 0]';
% x_new = x_old - inv(Df(x_old))*f(x_old);
% while norm(x_new - x_old, Inf) > 1e-6
%     x_old = x_new;
%     x_new = x_old - inv(Df(x_old))*f(x_old);
% end

% y = @(x) (4 -tan(x))/(3*x*x  + sin(x));
% f = @(x) (x^3 - 2*x*y(x) + y(x)^7 - 4*y(x)*x^3 - 5);
% syms x;
% derf = matlabFunction(diff(x^3 - 2*x*y(x) + y(x)^7 - 4*y(x)*x^3 - 5));
% init = 1;
% root = Newton(f,derf,init);

% f = @(x) (cos(x + sqrt(2)) + x*(x/2 + sqrt(2)));
% syms x;
% derf = matlabFunction(diff(df));
% init = 0;
% root = Newton(f,derf,init);
% mroot  = modifiedNewton(f,derf,init,4);

% f = @(x) (exp(6*x) + 3*exp(2*x)*(log(2))^2 - log(8)*exp(4*x) - (log(2))^3);
% syms x;
% derf = matlabFunction(diff(exp(6*x) + 3*exp(2*x)*(log(2))^2 - log(8)*exp(4*x) - (log(2))^3));
% init = 0;
% root = Newton(f,derf,init);
% mroot  = modifiedNewton(f,derf,init,3);
% 
% function root = Newton(f,derf,init)
%     i = 1; 
%     x(i) = init;
%     while abs(f(x(i))) > 1e-10
%         i = i + 1;
%         x(i) = x(i-1) - f(x(i-1))/derf(x(i-1));
%     end
%     root = x(i);
% end
% 
% function root  = modifiedNewton(f,derf,init,p)
%     i = 1;
%     x(i) = init;
%     while abs(f(x(i))) > 1e-10
%         i = i + 1;
%         x(i) = x(i-1) - (p*f(x(i-1)))/derf(x(i-1));
%     end
%     root = x(i);
% end

%Forward difference interpolation
f = @(x) exp(x);
x = [];
fx = [];
inc = 1/32;
x(1)=1;
fx(1) = f(x(1));
n = 101;

mu = @(x) (x-1)/(inc);

for i=2:n
    x(i) = x(i-1) + inc;
    fx(i) = f(x(i));
end

mat = zeros(n,n+1);

for i=1:n
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j =3:n+1
    for i = j-1:n
        mat(i,j) = mat(i,j-1)- mat(i-1,j-1);
    end
end

Forward_matrix = mat;
Forward = interpolate(2.25,mat,n,mu)
xx = [1.25:0.01:2.5];
yy = [];
i=1;
for x= [1.25:0.01:2.5]
    yy(i) = interpolate(x,mat,n,mu);
    i = i + 1;
end

figure(1);
plot(xx, exp(xx),'--' , 'LineWidth' ,2);
hold on;
plot(xx, yy);

%Backward difference interpolation

Backward_value = interpolate_backward(2.25,mat,n,inc)
xx = [1.25:0.01:2.5];
yy = [];
i=1;
for x= [1.25:0.01:2.5]
    yy(i) = interpolate_backward(x,mat,n,inc);
    i = i + 1;
end
figure(2);
plot(xx(75:end), exp(xx(75:end)),'--' , 'LineWidth' ,2);
hold on;
plot(xx(75:end), yy(75:end));

function y = interpolate_backward(x, mat, n, inc)
    y = mat(n,2);
    mu = (x-mat(n,1))/(inc);
%     for i = 1:n-1
%         req = mat(n,i+2);
%         bincoeff = nck(mu+i-1,i);
%         y = y + bincoeff*req;
%     end
    for i = 1:n-1
        f = mu;
        for j = 1:i-1
            f = f*(mu+j);
        end
        y = y + (f*mat(n,i+2)/factorial(i));
    end
end


function y = interpolate(x,mat,n,mu)
    y = 0;
    k = mu(x);
    for i = 1:n
        req = mat(i,i+1);
        bincoeff = nck(k,i-1);
        y = y + bincoeff*req;
    end
end

function ans = nck(x,i)
    res= [];
    for j = 0:i-1
        res = [res (x-j)];
    end
    ans = prod(res);
    ans = ans/factorial(i);
end

