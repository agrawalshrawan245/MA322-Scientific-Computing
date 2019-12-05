clear all;
clc;

% f = @(x) (x^5 - 7);
% 
% x(1) = 1;

%Part (i)
% for i = 1:4
%     x(i+1) = x(i)*(1 + ((7 - x(i)^5)/(x(i)^2))^3);
% end

%Part(ii)
% for i = 1:4
%     x(i+1) = x(i) - (x(i)^5 - 7)/(x(i)^2);
% end

%Part(c)
% for i = 1:4
%     x(i+1) = x(i) - (x(i)^5 - 7)/(5*x(i)^4);
% end

%Partd
% for i = 1:4
%     x(i+1) = x(i) - (x(i)^5 - 7)/(12);
% end

% f = @(x) (exp(x) + 2^(-x) + 2*cos(x) -6);
% syms x;
% derf = matlabFunction(diff(exp(x) + 2^(-x) + 2*cos(x) -6));
% root = Newton(f,derf,1,2);

% f = @(x) (exp(x) - 3*x*x);
% syms x;
% derf = matlabFunction(diff(exp(x) - 3*x*x));
% root = Newton(f,derf,0,1);

% f = @(x) (exp(x) - 1.5 - atan(x));
% syms x;
% derf = matlabFunction(diff(exp(x) - 1.5 - atan(x)));
% rootN = Newton(f,derf,-4,1);
% rootS = Secant(f,-10,2);

f = @(x) (x^3 - 12*x^2 + 3*x + 1);
syms x;
derf = matlabFunction(diff(x^3 - 12*x^2 + 3*x + 1));
rootN = Newton(f,derf,-4,1);
rootS = Secant(f,-1,2);

function root = Newton(f, derf, a, b)
    x(1) = (a+b)/2;
    i = 1;
    while abs(f(x(i))) > 1e-5
        i = i + 1;
        x(i) = x(i-1) - f(x(i-1))/derf(x(i-1));
    end
    root = x(i);
end

function root = Secant(f, a, b)
    x(1) = a;
    x(2) = b;
    i = 2;
    while abs(f(x(i))) > 1e-5
        i = i + 1;
        x(i) = x(i-1) - f(x(i-1))*((x(i-1) - x(i-2))/(f(x(i-1)) - f(x(i-2))));
    end
    root = x(i);
end
