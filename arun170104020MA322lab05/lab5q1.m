clear all;
clc;

syms x;
F = log(exp(x) +2);
df = matlabFunction(diff(F));

f = @(x) log(exp(x) +2);

n=4;
x = [-1 -1 -0.5 -0.5 0 0 0.5 0.5];
fx = [0.86199480 0.86199480 0.95802009 0.95802009 1.0986123 1.0986123 1.2943767 1.2943767];
dfx = [0.15536240 0.15536240 0.23269654 0.23269654 0.33333333 0.33333333 0.45186776 0.45186776];

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

val = 0.25;
pre = [(val-x(1))];
for i=2:n-1
    pre(i) = pre(i-1)*(val-x(i));
end

ans = fx(1);
for i =1:n-1
   ans = ans + pre(i)*mat(i+1,i+2);
end

err = abs(f(val) - ans);
fprintf('Value of f(0.25) is %.5f.\n',ans);
fprintf('The absolute error is %d.\n',err);
