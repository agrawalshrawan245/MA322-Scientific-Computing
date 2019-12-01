clear all;
clc;
clf;
f = @(x,y) (x+y);
sol = @(x) 2*exp(x) - x - 1;

prev = 1;
h = 0.005;
n = 1/0.005;


%By Euler's Method
fprintf('By Euler method:\n'); 
fprintf('t \t\t Approx \t Exact \t\t Error\n');
t = 0;
for i = 1:n
    y(i) = prev + h*f(t,prev);
    prev = y(i);
    t = t + h;
    err(i) = abs(sol(t) - y(i));
    fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));
end
figure(1);
plot(0.005:0.005:1, y);
hold on;
xx = 0.005:0.005:1;
plot(xx,sol(xx),'--', 'linewidth', 2);
legend('Approximation', 'Exact Solution');
xlabel('t(i)');
ylabel('y(t(i))')
figure(2);
plot(0.005:0.005:1, err);
xlabel('t(i)');
ylabel('Error');
figure(3);
semilogx(1:n,err);
title('Semilogx plot of t(i) vs Error');

%By Modified Euler's method
fprintf('\nBy modified Euler method:\n'); 
fprintf('t \t\t Approx \t Exact \t\t Error\n');
t = 0;
prev = 1;
for i = 1:n
    y(i) = prev + (h/2)*(f(t,prev) + f(t+h,prev+h*(f(t,prev))));
    prev = y(i);
    t = t + h;
    err(i) = abs(sol(t) - y(i));
    fprintf('%f \t %f \t %f \t %f\n', t,y(i),sol(t),err(i));
end
figure(4);
plot(0.005:0.005:1, y);
hold on;
xx = 0.005:0.005:1;
plot(xx,sol(xx),'--', 'linewidth', 2);
legend('Approximation', 'Exact Solution');
xlabel('t(i)');
ylabel('y(t(i))')
figure(5);
plot(0.005:0.005:1, err);
xlabel('t(i)');
ylabel('Error');
figure(6);
semilogx(1:n,err);
title('Semilogx plot of t(i) vs Error');
