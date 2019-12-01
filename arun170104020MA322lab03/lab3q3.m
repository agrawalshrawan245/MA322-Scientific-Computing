g = @(x) (4-tan(x))/(3*x^2 + sin(x));
f = @(x) (x^3 - 2*x*g(x) + (g(x))^7 - 4*(x^3)*g(x) - 5);
f2 = @(x) g(x)*sin(x) + 3*x*x*g(x) + tan(x) - 4;

syms x;
F = (x^3 - 2*x*g(x) + (g(x))^7 - 4*(x^3)*g(x) - 5);
df = diff(F);
derf = matlabFunction(df);

x = 1;
i =0;
while(abs(f(x)) > 1e-7)
    X(i+1,1) = i+1;
    X(i+1,2) = x;
    X(i+1,3) = g(x);
    X(i+1,4) = f(x);
    X(i+1,5) = f2(x);
    x = x - f(x)/(derf(x));
    i = i+1;
end
 X(i+1,1) = i;
 X(i+1,2) = x;
 X(i+1,3) = g(x);
 X(i+1,4) = f(x);
 X(i+1,5) = f2(x);