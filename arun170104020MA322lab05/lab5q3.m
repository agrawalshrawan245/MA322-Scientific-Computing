clear all;
clc;
n=5;
x = [0 0 3 3 5 5 8 8 13 13];
fx = [0 0 225 225 383 383 623 623 993 993];
dfx = [75 75 77 77 80 80 74 74 72 72];

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

syms val
pre = [(val-x(1))];
for i=2:n-1
    pre(i) = pre(i-1)*(val-x(i));
end

pos = fx(1);
for i =1:n-1
   pos = pos + pre(i)*mat(i+1,i+2);
end

speed = diff(pos);
acc = diff(speed);
acc = matlabFunction(acc);
speed = matlabFunction(speed);
pos = matlabFunction(pos);

fprintf('By Hermite interpolation , position of the car at t = 10 seconds is %.5f feet.\n', pos(10)); 
fprintf('By Hermite interpolation , speed of the car at t = 10 seconds is %.5f feet/sec.\n', speed(10));

figure(1);
xx =  [0:0.01:13];
plot(xx,speed(xx), 'g');
xlabel('Time in seconds');
ylabel('Speed in ft/s');
title('Plot of speed of the car vs. time');

fprintf('We have 1 mile/hr = 1.46666667 ft/s.\n');
target_speed = 55*1.46666667;
fprintf('Therefore 55 mi/h = %.5f ft/s.\n', target_speed);
 
first_time = 0;
max_speed = 0;
for t = [0:0.01:13]
    if speed(t)>=target_speed
        first_time = t;
        break;
    end
end

for t = [0:0.001:13]
    if speed(t)>max_speed
        max_speed = speed(t);
    end
end

fprintf('First time the car exceeds 55mi/h is %.2f s.\n', first_time);
fprintf('Predicted maximum speed of the car is %.3f ft/s.\n', max_speed);
