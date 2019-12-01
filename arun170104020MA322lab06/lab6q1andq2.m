clear all;
clc;


%Part a

syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
df = matlabFunction(df);
f = @(x) (2.*x)./(x.^2 - 4);
a = 1;
b = 1.6;

range = a:0.01:b;
mxdf = max(abs(df(range)));

I  = integral(f,a,b);
val = (b-a)*f(a);
fprintf('By Rectangle rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf*(b-a)*(b-a))/2;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

%Part b
syms x;
df = exp(3*x)*sin(2*x);
df = diff(df);
df = matlabFunction(df);
f = @(x) exp(3.*x).*sin(2.*x);
a = 0;
b = pi/4;

range = a:0.01:b;
mxdf = max(abs(df(range)));

I = integral(f,a,b);
val = (b-a)*f(b);
fprintf('By Rectangle rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf*(b-a)*(b-a))/2;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

%Part c
syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
df = matlabFunction(df);
f = @(x) (sin(x).*sin(x) - 2.*x.*sin(x) + 1);
a = 0.75;
b = 1.3;

range = a:0.01:b;
mxdf = max(abs(df(range)));

I = integral(f,a,b);
val = (b-a)*f(a);
fprintf('By Rectangle rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf*(b-a)*(b-a))/2;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

%Part d
syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
df = matlabFunction(df);
f = @(x) 1./(x.*log(x));
a = exp(1);
b = exp(1) + 1;

range = a:0.01:b;
mxdf = max(abs(df(range)));

I = integral(f,a,b);
val = (b-a)*f(a);
fprintf('By Rectangle rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf*(b-a)*(b-a))/2;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);
