clear all;
close all;

f1 = @(x) 3*x.^2 - exp(x);
f2 = @(x) x - cos(x);
g1 = @(x) log(3*x.^2);
g2 = @(x) cos(x)

x(1)=1;
i=1;
xp = [x(1)];
yp = [x(1)];
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f1(x(1)));
while(abs(f1(x(i)))> 1e-3)
    i =i+1;
    x(i) = g1(x(i-1));
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f1(x(i)));
    xp=[xp x(i-1)];
    xp=[xp x(i)];
    yp=[yp x(i)];
    yp=[yp x(i)];
end

root = x(i)

hold on;
fplot(g1,[1 4]);
fplot(@(x) x,[1 4]);
plot(xp,yp,'k');
hold off

figure;
x(1)=1;
i=1;
xp = [x(1)];
yp = [x(1)];
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f2(x(1)));
while(abs(f2(x(i)))> 1e-3)
    i =i+1;
    x(i) = cos(x(i-1));
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f2(x(i)));
    xp=[xp x(i-1)];
    xp=[xp x(i)];
    yp=[yp x(i)];
    yp=[yp x(i)];
end

root = x(i)

hold on;
fplot(g2,[-0.5 1.5]);
fplot(@(x) x,[-0.5 1.5]);
plot(xp,yp,'k');
hold off

