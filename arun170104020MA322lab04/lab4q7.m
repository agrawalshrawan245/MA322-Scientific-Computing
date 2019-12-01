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