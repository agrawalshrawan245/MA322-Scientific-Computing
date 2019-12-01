clear all;
clc;

syms x;
df = 4/(1 + x^2);
df = diff(df);
df = matlabFunction(df);

f = @(x) 4/(1 + x^2);

a = 0;
b = 1;
h = b-a;
val1 = f(a)*(b-a);
fprintf('Using Rectangle rule, the estimated value is %.5f and absolute error is %.5f.\n', val1, abs(pi-val1));

val2 = ((f(a) + f(b))*h)/2;
fprintf('Using Trapezoidal rule, the estimated value is %.5f and absolute error is %.5f.\n', val2, abs(pi-val2));

val3 = val2 - ((h*h)/12)*(df(b) - df(a));
fprintf('Using Corrected Trapezoidal rule, the estimated value is %.5f and absolute error is %.5f.\n', val3, abs(pi-val3));

val4 = ((b-a)/6)*(f(a) + 4*f((a + b)/2) + f(b));
fprintf('Using Simpsons one-third rule, the estimated value is %.5f and absolute error is %.5f.\n', val4, abs(pi-val4));

inc = (b-a)/3;
x = [a a+inc a+2*inc a+3*inc];
val5 = ((3*inc)/8)*(f(x(1)) + 3*f(x(2)) + 3*f(x(3)) + f(x(4)));
fprintf('Using Simpsons three-eighth rule, the estimated value is %.5f and absolute error is %.5f.\n', val5, abs(pi-val5));
