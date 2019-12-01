clear all;
clc;
clf;

% Given PDE as ut + a*ux = 0 . We have to determine solution for three time
% levels

a = -2;
initial_cond = @(x) (1 + sin(2*pi*x));
xl = 0;
xr = 1;
h = 0.005;
k = 0.001;
Courant_no = (a*k)/h;
m = (xr - xl)/h;
u(1,:) = initial_cond(xl:h:xr);
u(1,m+1) = 1 ; %Right B.C.
for i = 1:m-1
    A(i,i) = 1 + Courant_no;
    A(i,i+1) = -Courant_no;
end
A(m,m) = 1 + Courant_no;

figure(1);
plot(xl:h:xr, u(1,:));
xlabel('x(i)');
ylabel('u(i)');
title('Numerical Solution at different time levels by FTFS');    
for i = 2:4
    u(i,(1:m)) = (A*(u(i-1,(1:m))'))';
    u(i,m) = u(i,m) + Courant_no*u(i-1,m+1);
    u(i,m+1) = 1;
    hold on;
    plot(xl:h:xr, u(i,:));
end

% By BTFS scheme
clear u A;
for i = 1:m
    A(i,i) = 1 - Courant_no;
    A(i,i+1) = Courant_no;
end
A(m+1,m+1) = 1;
u(1,:) = initial_cond(xl:h:xr);
figure(2);
plot(xl:h:xr, u(1,:));
xlabel('x(i)');
ylabel('u(i)');
title('Numerical Solution at different time levels by BTFS');
for i = 2:4
    u(i,:) = (A\(u(i-1,:)'))';
    hold on;
    plot(xl:h:xr, u(i,:));
end