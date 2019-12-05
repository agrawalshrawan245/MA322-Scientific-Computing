clear all;
clc;

f = @(x) (exp(x) - sin(x));
a = -4;
b = -2;
root = bisection(f, a, b);


f = @(x) (2 + cos(exp(x) - 2) - exp(x));
a = 0.5;
b = 1.5;
aroot = bisection(f, a, b);

f = @(x) (x - 2^(-x));
a = 0;
b = 1;
broot = bisection(f, a, b);


f = @(x) (exp(x) -x^2  + 3*x - 2);
a = 0;
b = 1;
croot = bisection(f, a, b);

f = @(x) (2*x*cos(2*x) - (1+x)^2);
droot1 = bisection(f, -3, -2);
droot2 = bisection(f, -1, 0);

f = @(x) (x*cos(x) - 2*x*x + 3*x - 1);
eroot1 = bisection(f, 0.2, 0.3);
eroot2 = bisection(f, 1.2, 1.3);

f = @(x) (x^3 - 25);
cuberoot = bisection(f, 2, 3);

function root = bisection(f, a, b)
    c = (a+b)/2;
    while abs(f(c))> 1e-4
      if (f(a)*f(c)) < 0
          b = c;
      else
          a = c;
      end
      c = (a+b)/2;
    end
    root = c;
end