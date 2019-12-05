clear all;
clc;

% Pde is of the form ut = alpha * uxx;

% exact_sol = @(x,t) (exp(-t).*sin((pi*x)./2) + exp(-t/4).*sin((pi*x)./4));
% initial_cond = @(x) sin((pi*x)./4).*(1 + 2*(cos((pi*x)./4)));
% alpha = (4/(pi*pi));
% bcl = @(t) 0;
% bcr = @(t) 0;
% 
% ti = 0;
% tf = 1;
% xi = 0;
% xf = 4;

exact_sol = @(x,t) (exp((-pi*pi*t)/4).*sin((pi*x)/2) + (1/2)*exp(-4*pi*pi*t).*sin(2*pi*x));
initial_cond = @(x) sin((pi*x)/2) + (1/2)*sin(2*pi*x);
alpha = 1;
bcl = @(t) 0;
bcr = @(t) exp((-pi*pi*t)/4);

ti = 0;
tf = 1;
xi = 0;
xf = 1;

h = 0.01;
k = 0.01;

x = xi:h:xf;
t = ti:k:tf;
n = length(t);
m = length(x);
% lambda = alpha*(k/(h^2));
% u(:,1) = initial_cond(x);

%FTCS
% A(1,1) = 1;
% for i = 2:m-1
%     A(i,i-1) = lambda;
%     A(i,i) = 1 - 2*lambda;
%     A(i,i+1) = lambda;
% end
% A(m,m) = 1;
% for i = 2:n
%     u(:,i) = A*(u(:,i-1));
%     u(1,i) = bcl(t(i));
%     u(m,i) = bcr(t(i));
% end
% 
% figure(1);
% [X,T] = meshgrid(x,t);
% Z = exact_sol(X,T);
% surf(X,T,Z);
% figure(2);
% surf(X,T,u');
% figure(3);
% plot(x,u(:,n)','--','LineWidth',2);
% hold on;
% plot(x, exact_sol(x,1));

%BTCS
% A(1,1) = 1;
% for i = 2:m-1
%     A(i,i-1) = -lambda;
%     A(i,i) = 1 + 2*lambda;
%     A(i,i+1) = -lambda;
% end
% A(m,m) = 1;
% for i = 2:n
%     tmp = u(:,i-1);
%     tmp(1,1) = bcl(t(i));
%     tmp(m,1) = bcr(t(i));
%     u(:,i) = A\tmp;
% end
% 
% figure(1);
% [X,T] = meshgrid(x,t);
% Z = exact_sol(X,T);
% surf(X,T,Z);
% figure(2);
% surf(X,T,u');
% figure(3);
% plot(x,u(:,n)','--','LineWidth',2);
% hold on;
% plot(x, exact_sol(x,tf));

%Crank Nicolson
% lambda = alpha*(k/(h^2));
% u(:,1) = initial_cond(x);
% A1(1,1) = 1;
% A1(m,m) = 1;
% for i = 2:m-1
%     A1(i,i-1) = lambda/2;
%     A1(i,i) = (1 - lambda);
%     A1(i,i+1) = lambda/2;
% end
% 
% A0(1,1) = 1;
% A0(m,m) = 1;
% for i = 2:m-1
%     A0(i,i-1) = -lambda/2;
%     A0(i,i) =  1 + lambda;
%     A0(i,i+1) = -lambda/2;
% end
% 
% for i = 2:n
%     b = A1*u(:,i-1);
%     b(1,1) = bcl(t(i));
%     b(m,1) = bcr(t(i));
%     u(:,i) = A0\b;
% end
% 
% figure(1);
% [X,T] = meshgrid(x,t);
% Z = exact_sol(X,T);
% surf(X,T,Z);
% figure(2);
% surf(X,T,u');
% figure(3);
% plot(x,exact_sol(x,tf));
% hold on;
% plot(x, u(:,n)', '--', 'LineWidth',2);


