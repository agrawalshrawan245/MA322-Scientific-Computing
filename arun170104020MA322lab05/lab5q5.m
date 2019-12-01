clear all;
clc;

f = @(x) x*cos(x) - 2*x*x + 3*x - 1;

inc= 0.1;
a = 0.1;
b = 0.4;
x = [a:inc:b];
n = length(x);
y = [-0.62049958 -0.28398668 0.00660095 0.24842440];

for i=1:n-1
    h(i) = inc;
end

dya = 3.58502082;
dyb = 2.16529366;

mat = zeros(n,n);

mat(1,1) = h(1)/3;
mat(1,2) = h(1)/6;
for i=2:n-1
    mat(i,i-1) = h(i-1)/6;
    mat(i,i) = (h(i-1) + h(i))/3;
    mat(i,i+1) = h(i)/6;
end 
mat(n,n-1) = h(n-1)/6;
mat(n,n) = h(n-1)/3;

d(1) = ((y(2)-y(1))/h(1))-dya;
for i=2:n-1
    d(i) = ((y(i+1)-y(i))/h(i)) -((y(i)-y(i-1))/h(i-1));
end
d(n) = dyb - ((y(n)-y(n-1))/h(n-1));
d =  d';

M = inv(mat)*d;
x_val = 0.2013;

ans = (M(2)*(x(3)-x_val)^3 + M(3)*(x_val-x(2))^3)/(6*h(2));
ans = ans + (x(3)-x_val)*((y(2)/h(2))-((h(2)*M(2))/6));
ans = ans + (x_val-x(2))*((y(3)/h(2))-((h(2)*M(3))/6));

Exact = f(x_val);
Ans = ans;

Error = abs(Exact - Ans);

fprintf('Value of f(0.2013) is %.5f.\n',Ans);
fprintf('The absolute error is %d.\n',Error);