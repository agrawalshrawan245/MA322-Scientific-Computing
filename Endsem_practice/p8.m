clear all;
clc;

%given pde as ut + aux = 0;
% initial_cond = @(x) 1 + sin(2*pi*x);
% bcr = @(t) 1;
% a = -2;
% xi = 0;
% xf = 1;
% h = 0.005;
% k = 0.001;
% x = xi:h:xf;
% m = length(x);
% u(:,1) = initial_cond(x);
% R = a*(k/h);
% %FTFS
% for i = 1:m-1
%     A(i,i) = 1 + R;
%     A(i,i+1) = -R;
% end
% A(m,m) = 1 + R;
% plot(x,u(:,1)');
% for i = 2:4
%     u(:,i) = A*u(:,i-1);
%     u(m,i) = u(m,i) - R*bcr((i-2)*k);
%     hold on;
%     plot(x,u(:,i)');
% end

%BTFS
% for i =1:m-1
%     A(i,i) = 1 - R;
%     A(i,i+1) = R;
% end
% A(m,m) = 1- R;
% plot(x,u(:,1)');
% for i = 2:4
%     b = u(:,i-1);
%     b(m,1) = b(m,1) - R*bcr((i-1)*h);
%     u(:,i) = A\b;
%     hold on;
%     plot(x,u(:,i)');
% end

% initial_cond = @(x) (1 + sin(2*pi*x));
% xi = 0;
% xf = 1;
% h = 0.005;
% k = 0.001;
% x = xi:h:xf;
% m = length(x);
% a = 2;
% R = a*(k/h);
% bcl = @(t) 1;
% 
% u(:,1) = initial_cond(x);
% plot(x,u(:,1)');
% A(1,1) = 1;
% A(m,m) = 1;
% for i = 2:m-1
%     A(i,i-1) = (R^2 + R)/2;
%     A(i,i) = 1 - R^2;
%     A(i,i+1) = (R^2 - R)/2;
% end
% 
% for i = 2:4
%     b = u(:,i-1);
%     u(:,i) = A*b;
%     u(1,i) = bcl((i-1)*k);
%     u(m,i) = bcl((i-1)*k);
%     hold on;
%     plot(x,u(:,i)');
% end

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









clear all;
clc;
x = [2 3 5 6];
fx = [1.5713 1.5719 1.5738 1.5751];

size = 4;
mat = zeros(size,size+1);

for i = 1:size
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j = 3:size+1
    for i = j-1:size
        mat(i,j) = (mat(i,j-1) - mat(i-1,j-1))/(mat(i,1) - mat(i-j+2,1));
    end
end


%Calculation of f(4)
x_val = 4;
pre = [(x_val-x(1))];
for i=2:size-1
    pre(i) = pre(i-1)*(x_val-x(i));
end

%Using second degree interpolating polynomial
y_val1 = fx(1);
for i =1:size-2
    y_val1 = y_val1 + pre(i)*mat(i+1,i+2);
end

y_val1

%Using third degree interpolating polynomial
y_val2 = fx(1);
for i =1:size-1
    y_val2 = y_val2 + pre(i)*mat(i+1,i+2);
end

y_val2