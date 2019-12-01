clear all;
clc;

N = 4;
a=1;
b=2;
h = (b-a)/N;
x = [a:h:b];
y = [10 8 7 6 5];

%Using formula for composite trapezoidal rule

val = y(1) + y(N+1);

for i=2:N
 val = val + 2*y(i);  
end

val  = val*(h/2);

fprintf('Estimate by the composite trapezoid rule is: %0.5f.\n', val);
