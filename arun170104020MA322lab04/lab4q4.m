%Lagrange polynomial interpolation
n=5;
x = [0 0.1 0.3 0.6 1.0];
fx = [-6 -5.89483 -5.65014 -5.17788 -4.28172];

ans =0;

for i=1:n
    tmp = 1;
    for j =1:n
        if i~=j       
            tmp = tmp*( (0.2-x(j))/(x(i)-x(j)) );
        end
    end
    ans = ans + tmp*fx(i);
end

fprintf('f(0.2) by Lagrange interpolation= %f\n', ans);

x = [x 1.1];
fx = [fx -3.99583];
n=6;

ans =0;

for i=1:n
    tmp = 1;
    for j =1:n
        if i~=j       
            tmp = tmp*( (0.2-x(j))/(x(i)-x(j)) );
        end
    end
    ans = ans + tmp*fx(i);
end

fprintf('f(0.2) by Lagrange interpolation adding_extra_point = %f\n', ans);

%Newton's divided difference interpolation

x = [0 0.1 0.3 0.6 1.0];
fx = [-6 -5.89483 -5.65014 -5.17788 -4.28172];

size = 5;
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

x_val = 0.2;
pre = [(x_val-x(1))];
for i=2:size-1
    pre(i) = pre(i-1)*(x_val-x(i));
end

y_val = fx(1);
for i =1:size-1
    y_val = y_val + pre(i)*mat(i+1,i+2);
end

fprintf('f(0.2) by Newton divided difference interpolation = %f\n', y_val);

x = [x 1.1];
fx = [fx -3.99583];
size=6;
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

x_val = 0.2;
pre = [(x_val-x(1))];
for i=2:size-1
    pre(i) = pre(i-1)*(x_val-x(i));
end

y_val = fx(1);
for i =1:size-1
    y_val = y_val + pre(i)*mat(i+1,i+2);
end

fprintf('f(0.2) by Newton divided difference interpolation adding extra point = %f\n', y_val);
