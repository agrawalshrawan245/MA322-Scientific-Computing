clear all;
clc;
clf;
f = @(x) atan(x);

n = 11;
a = 1;
b = 6;
inc = (b-a)/(n-1);
x = [a:inc:b];
fx = f(x);
mu = @(x) (x-a)/inc;

mat = zeros(n,n+1);

for i=1:n
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j = 3:n+1
    for i = j-1:n
        mat(i,j) = mat(i,j-1) - mat(i-1,j-1);
    end
end

syms y;
p = fx(1);

disp('Coefficients in the Newton form of polynomial are:');
for i=1:n-1
    p = p + nck(mu(y), i)*mat(i+1,i+2);
    fprintf('%0.5f\n', mat(i+1,i+2));
end

p = matlabFunction(p);

%Computing and printing the difference between the polynomial and the
%function at specified points
xx = [0:0.25:8];
yy = f(xx);
 

plot(xx,yy,'LineWidth', 2);
hold on;
plot(xx,p(xx), 'g');
title('Green curve is interpolant and blue is arctan(x)');


diff = [p(xx).' f(xx).' (p(xx) - f(xx)).'];
fprintf('The difference at required points is in table below:\n');
fprintf('\tp(x) \t\t f(x) \t\t Difference(p(x) - f(x))\n');

disp(compose('%.10f',diff));

function out = nck(n,k)
    tmp = 1;
    for i = 0:k-1
        tmp = tmp * (n-i);
    end
    tmp = tmp/factorial(k);
    out = tmp;
end