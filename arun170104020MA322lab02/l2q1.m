close all;
clear all;
f = @(x) x.^3 - x - 1;
g = @(x) (x+1).^(1/3);

x(1)=1;
i = 1;

xp=[x(1)];
yp=[x(1)];
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
while(abs(f(x(i)))> 1e-2)
    i = i + 1;
    x(i) = g(x(i-1));
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
    xp=[xp x(i-1)];
    xp=[xp x(i)];
    yp=[yp x(i)];
    yp=[yp x(i)];
end

root = x(i)

hold on;
fplot(g,[1 2]);
fplot(@(x) x,[1 2]);
plot(xp,yp,'k');
hold off




