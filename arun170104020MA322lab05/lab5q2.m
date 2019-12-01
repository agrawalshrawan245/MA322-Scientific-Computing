clear all;
clc;

syms x;
F = x*cos(x) - 2*x^2 + 3*x - 1;
df = matlabFunction(diff(F));

f = @(x) x*cos(x) - 2*x^2 + 3*x - 1;

n=4;
x = [0.1 0.1 0.2 0.2 0.3 0.3 0.4 0.4];
fx = [-0.62049958 -0.62049958 -0.28398668 -0.28398668 0.00660095 0.00660095 0.24842440 0.24842440];
dfx = [3.58502082 3.58502082 3.14033271 3.14033271 2.66668043 2.66668043 2.16529366 2.16529366];

n = 2*n;
mat = zeros(n, n+1);

for i=1:n
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j = 3:n+1
    for i = j-1:n
        if(mat(i,1)~= mat(i-j+2,1))
            mat(i,j) = (mat(i,j-1) - mat(i-1,j-1))/(mat(i,1) - mat(i-j+2,1));
        else
            mat(i,j) = dfx(i);
        end
    end
end

val = 0.2013;
pre = [(val-x(1))];
for i=2:n-1
    pre(i) = pre(i-1)*(val-x(i));
end

ans = fx(1);
for i =1:n-1
   ans = ans + pre(i)*mat(i+1,i+2);
end

err = abs(f(val) - ans);
fprintf('Value of f(0.2013) is %.5f.\n',ans);
fprintf('The absolute error is %d.\n',err);