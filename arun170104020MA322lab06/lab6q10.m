clear all;
clc;

f = @(x) 4*(((3.*cos(x)).^2 + (2.*sin(x)).^2).^(0.5));

a = 0;
b = pi/2;

h = (b - a)/2;
val = (h/3)*(f(a) + 4*f((a + b)/2) + f(b));

I = integral(f,0,pi/2,'AbsTol',1e-7);
while (abs(I-val) > 1e-6)
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

fprintf('The length of the graph of the ellipse is %.5f.\n', val);

