clear all;
clc;

actual_val = log(4) - (3/4);

f = @(x) x.*log(x);

a = 1;
b = 2;

%By composite trapezoidal rule
h = b-a;
val = (h/2)*(f(a) + f(b));
n=1;
while abs(actual_val - val) > 1e-5
    n = n + 1;
    h = (b-a)/n;
    x = [];
    val  = f(a) + f(b);
    for i =1:n-1
        x(i) = a + i*h;
        val = val + 2*f(x(i));
    end
    val = val*(h/2);
end
fprintf('By composite trapezoidal rule, the approximate value is %.6f.\n', val);
fprintf('Value of n is %d and value of h is %.5f.\n', n, h);

%By composite Simpson's rule
h = (b - a)/2;
val = (h/3)*(f(a) + 4*f((a + b)/2) + f(b));
while abs(actual_val-val) > 1e-5
    h = h/2;
    x = [a:2*h:b];
    n = length(x)-1;
    y = f(x);
    val = y(1) + y(n+1);
    for i=2:n
        val = val + 2*y(i);
    end
    for i=1:n
        val = val + 4*f((x(i) + x(i+1))/2);
    end
    val = val*(h/3);
end

fprintf('By composite simpsons rule, the approximate value is %.6f.\n', val);
fprintf('Value of n is %d and value of h is %.5f.\n', 2*n, h);

%By Composite midpoint rule
h = b-a;
val = f((a+b)/2)*h;
n=1;
while abs(actual_val-val) > 1e-5
    n = n+1;
    h = (b-a)/n;
    val = 0;
    for i=1:n
        val = val + f(a+(i-(1/2))*h);
    end
    val = val*h;
end

fprintf('By composite midpoint rule, the approximate value is %.6f.\n', val);
fprintf('Value of n is %d and value of h is %.5f.\n', 2*n, h);