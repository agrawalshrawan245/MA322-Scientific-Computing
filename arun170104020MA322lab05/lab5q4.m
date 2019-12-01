clear all;
clc;

x = [0.3:0.02:0.44];
y = sin(x);

n = length(x);

for i=1:n-1
    h(i) = 0.02;
end

mat = zeros(n-2,n-2);

mat(1,1) = (h(1) + h(2))/3;
mat(1,2) = h(2)/6;

for i=2:n-3
    mat(i,i-1) = h(i)/6;
    mat(i,i) = (h(i) + h(i+1))/3;
    mat(i,i+1) = h(i+1)/6;
end

mat(n-2,n-3) = h(n-2)/6;
mat(n-2,n-2) = (h(n-2) + h(n-1))/3;

for i=1:n-2
    d(i) = ((y(i+2)-y(i+1))/h(i+1)) - ((y(i+1) - y(i))/h(i));
end
d = d';
M = inv(mat)*d;

x_val = 0.3102;

ans = (((x_val - x(1))^3)*(M(1)))/(6*h(1));
ans = ans + (y(1)/h(1))*(x(2)-x_val);
ans = ans + ((y(2)/h(1)) - (h(1)*M(1))/6)*(x_val-x(1));


Exact = sin(x_val);
Ans = ans;

Error = abs(Exact - Ans);
fprintf('Value of f(0.3102) is %.5f.\n',Ans);
fprintf('The absolute error is %d.\n',Error);

