clear all;
clc;
clf;

f = @(x) 4./(1 + x.^2);
g = @(x) sqrt(1-x.^2) - x;

fa = 0;
fb = 1;

ga = 0;
gb = 1/sqrt(2);


%Code for function f
h = (fb - fa)/2;
val = (h/3)*(f(fa) + 4*f((fa + fb)/2) + f(fb));
subint = [2];
err = [abs(pi-val)];
ind = 1;
while abs(pi-val) > (1e-5)/2
    h = h/2;
    ind = ind + 1;
    subint(ind) = subint(ind-1)*2;
    x = [fa:2*h:fb];
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
    err(ind) = abs(pi-val);
end

fprintf('Approximate value using Simpsons one third rule is %.5f.\n', val);
figure(1);
semilogy(subint,err);
title('Sub-intervals vs. error plot');

%Code for function g
h = (gb - ga)/2;
val = (h/3)*(g(ga) + 4*g((ga + gb)/2) + g(gb));
subintg = [2];
gerr = [abs((pi/8)-val)];
ind = 1;
while abs((pi/8)-val) > (1e-5)/2
    h = h/2;
    ind = ind + 1;
    subintg(ind) = subintg(ind-1)*2;
    x = [ga:2*h:gb];
    n = length(x)-1;
    y = g(x);
    val = y(1) + y(n+1);
    for i=2:n
        val = val + 2*y(i);
    end
    for i=1:n
        val = val + 4*g((x(i) + x(i+1))/2);
    end
    val = val*(h/3);
    gerr(ind) = abs((pi/8)-val);
end
fprintf('Approximate value using Simpsons one third rule is %.5f.\n', val);
figure(2);
semilogy(subintg,gerr);
title('Sub-intervals vs. error plot');

