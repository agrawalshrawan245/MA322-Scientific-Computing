clc;
close all;
clear all;

f1 = @(x) 6*x(1) - 2*cos(x(2)*x(3)) - 1;
f2 = @(x) 9*x(2) + sqrt((x(1)^2 + sin(x(3)) + 1.06)) + 0.9;
f3 = @(x) 60*x(3) + 3*exp(-x(1)*x(2)) + 10*pi - 3;

f = @(x) [6*x(1) - 2*cos(x(2)*x(3)) - 1; 9*x(2) + sqrt((x(1)^2 + sin(x(3)) + 1.06)) + 0.9; 60*x(3) + 3*exp(-x(1)*x(2)) + 10*pi - 3];

Df = @(x) [6, 2*x(3)*sin(x(2)*x(3)), 2*x(2)*sin(x(2)*x(3)); (x(1)/sqrt(x(1)^2 + sin(x(3)) + 1.06)), 9, (cos(x(3))/(2*sqrt(x(1)^2 + sin(x(3)) + 1.06))); -3*x(2)*exp(-x(1)*x(2)), -3*x(1)*exp(-x(1)*x(2)), 60];

xold = [0;0;0];

x = xold - inv(Df(xold))*f(xold) ;
X(1,1) = 1;
X(1,2) = x(1);
X(1,3) = x(2);
X(1,4) = x(3);
X(1,5) = f1(x);
X(1,6) = f2(x);
X(1,7) = f3(x);
i=1;
while(norm((x - xold), Inf) > 1e-6)
    tmp = x;
    x = x - inv(Df(x))*f(x) ;
    i = i+1;
    X(i,1) = i;
    X(i,2) = x(1);
    X(i,3) = x(2);
    X(i,4) = x(3);
    X(i,5) = f1(x);
    X(i,6) = f2(x);
    X(i,7) = f3(x);
    xold = tmp;
end
