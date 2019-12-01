clear all;
clc;

x = [0 3 5 8 13];
y = [0 225 383 623 993];

n = length(x);

for i=1:n-1
    h(i) = x(i+1)-x(i);
end

%Natural Spline interpolation
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

x_val = 10;

pos = (M(3)*(x(5)-x_val)^3)/(6*h(4));
pos = pos + ((y(4)/h(4)) - (h(4)*M(3))/6)*(x(5)-x_val);
pos = pos + (y(5)/h(4))*(x_val-x(4));

Position_at_10 = pos;

speed = -(M(3)*(x(5)-x_val)^2)/(2*h(4));
speed = speed + (y(5)-y(4))/h(4);
speed = speed + (M(3)*h(4))/6;

Speed_at_10 = speed;

fprintf('By natural spine interpolation , position of the car at t = 10 seconds is %.5f feet.\n', Position_at_10); 
fprintf('By natural spine interpolation , speed of the car at t = 10 seconds is %.5f feet/sec.\n', Speed_at_10); 

%Clamped Spline interpolation

dyi = 75;
dyf = 72;

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

d(1) = ((y(2)-y(1))/h(1)) - dyi;
for i=2:n-1
    d(i) = ((y(i+1)-y(i))/h(i)) -((y(i)-y(i-1))/h(i-1));
end
d(n) = dyf - ((y(n)-y(n-1))/h(n-1));

M = inv(mat)*d;
x_val = 10;

pos = (M(4)*(x(5)-x_val)^3 + M(5)*(x_val-x(4))^3)/(6*h(4));
pos = pos + (x(5)-x_val)*((y(4)/h(4))-((h(4)*M(4))/6));
pos = pos + (x_val-x(4))*((y(5)/h(4))-((h(4)*M(5))/6));

Clamped_Position_at_10 = pos;

speed = (-(M(4)*(x(5)-x_val)^2) + (M(5)*(x_val-x(4))^2))/(2*h(4));
speed = speed + (y(5)-y(4))/h(4);
speed = speed - ((M(5)-M(4))*h(4))/6;

Clamped_Speed_at_10 = speed;

fprintf('By clamped spine interpolation , position of the car at t = 10 seconds is %.5f feet.\n',Clamped_Position_at_10); 
fprintf('By clamped spine interpolation , speed of the car at t = 10 seconds is %.5f feet/sec.\n', Clamped_Speed_at_10); 
