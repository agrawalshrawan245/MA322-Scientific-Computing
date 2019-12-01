clear all;
clc;
clf;
n = 4;
x = [1 1.1 1.3 1.4];
fx =[];
for i = 1:n
    fx(i) = log(x(i));
end

syms y;
f = 0;
for i = 1:n
    num = 1;
    denom =1;
    for j = 1:n
      if i~=j
       num = num *(y-x(j));
       denom  = denom*(x(i) - x(j));
      end
    end
    f = f + fx(i)*(num/denom);
end

f =  matlabFunction(f);

figure(1);
plot(x, f(x), '--', 'LineWidth', 2);
hold on;
plot(x,log(x), 'g', 'LineWidth', 1); 
title('Blue -- is interpolant and ln(x) is in green colour');

max = -inf;
for x = [1:0.01:1.4]
    if abs(f(x) - log(x)) > max
        max = abs(f(x) - log(x));
    end
end

bound_for_absolute_error = max
figure(2);
err = [];
for x = [1:0.01:1.4]
    err = [err abs(f(x) - log(x))];
end

loglog([1:0.01:1.4], err);

