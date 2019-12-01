clear all;
clc;
clf;
%Forward difference interpolation
f = @(x) exp(x);
x = [];
fx = [];
inc = 1/32;
x(1)=1;
fx(1) = f(x(1));
n = 101;

mu = @(x) (x-1)/(inc);

for i=2:n
    x(i) = x(i-1) + inc;
    fx(i) = f(x(i));
end

mat = zeros(n,n+1);

for i=1:n
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j =3:n+1
    for i = j-1:n
        mat(i,j) = mat(i,j-1)- mat(i-1,j-1);
    end
end

Forward_matrix = mat;
Forward = interpolate(2.25,mat,n,mu)
xx = [1.25:0.01:2.5];
yy = [];
i=1;
for x= [1.25:0.01:2.5]
    yy(i) = interpolate(x,mat,n,mu);
    i = i + 1;
end

figure(1);
plot(xx, exp(xx),'--' , 'LineWidth' ,2);
hold on;
plot(xx, yy);

%Backward difference interpolation

Backward_value = interpolate_backward(2.25,mat,n,inc)
xx = [1.25:0.01:2.5];
yy = [];
i=1;
for x= [1.25:0.01:2.5]
    yy(i) = interpolate_backward(x,mat,n,inc);
    i = i + 1;
end
figure(2);
plot(xx, exp(xx),'--' , 'LineWidth' ,2);
hold on;
plot(xx, yy);

function y = interpolate_backward(x, mat, n, inc)
    y = mat(n,2);
    mu = (x-mat(n,1))/(inc);
%     for i = 1:n-1
%         req = mat(n,i+2);
%         bincoeff = nck(mu+i-1,i);
%         y = y + bincoeff*req;
%     end
    for i = 1:n-1
        f = mu;
        for j = 1:i-1
            f = f*(mu+j);
        end
        y = y + (f*mat(n,i+2)/factorial(i));
    end
end


function y = interpolate(x,mat,n,mu)
    y = 0;
    k = mu(x);
    for i = 1:n
        req = mat(i,i+1);
        bincoeff = nck(k,i-1);
        y = y + bincoeff*req;
    end
end

function ans = nck(x,i)
    res= [];
    for j = 0:i-1
        res = [res (x-j)];
    end
    ans = prod(res);
    ans = ans/factorial(i);
end

