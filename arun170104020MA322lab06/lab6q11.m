clear all;
clc;
clf;

t = [0:6:84];
s = [124 134 148 156 147 133 121 109 99 85 78 89 104 116 123];

n = length(t);

syms y;
f = 0;
for i = 1:n
    num = 1;
    denom =1;
    for j = 1:n
      if i~=j
       num = num *(y-t(j));
       denom  = denom*(t(i) - t(j));
      end
    end
    f = f + s(i)*(num/denom);
end

f =  matlabFunction(f);

xx = [0:0.01:84];
plot(xx, f(xx));
%Using composite simpson's rule to calculate the length of the track
h = 6;
x = [t(1):2*h:t(n)];
n1 = length(x)-1;
y = f(x);
val = y(1) + y(n1+1);
for i=2:n1
    val = val + 2*y(i);
end
for i=1:n1
    val = val + 4*f((x(i) + x(i+1))/2);
end
val = val*(h/3);
fprintf('Using composite simpsons rule for calculation,  the length of the track obtained is %.2f feet.\n', val);
