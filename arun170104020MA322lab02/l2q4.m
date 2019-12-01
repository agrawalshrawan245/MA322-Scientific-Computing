clear all;
close all;

g = @(x) x.^2 + x - 0.75;
f = @(x) x.^2 - 0.75;

x(1)  = -0.8;
i=1;
xp = [x(1)];
yp = [x(1)];
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
while(abs(f(x(i)))> 1e-3)
    i = i+1;
    x(i) = g(x(i-1));
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
    xp=[xp x(i-1)];
    xp=[xp x(i)];
    yp=[yp x(i)];
    yp=[yp x(i)];
end

hold on;
title('Initial value -0.8');
fplot(g,[-1 -0.5]);
fplot(@(x) x,[-1 -0.5]);
plot(xp,yp,'k');
hold off

disp('With initial value -0.8 solution is = ' + string(x(i)));

y(1)  = 0.8;
i=1;
xp = [y(1)];
yp = [y(1)];
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(1)));
while(abs(f(y(i)))> 1e-3)
    i = i+1;
    y(i) = g(y(i-1));
    fprintf('\t %d \t\t %f \t %f\n',i,y(i),f(y(i)));
    xp=[xp y(i-1)];
    xp=[xp y(i)];
    yp=[yp y(i)];
    yp=[yp y(i)];
end

disp('With initial value 0.8 solution is = ' + string(y(i)));

figure;
hold on;
title('Initial value 0.8');
fplot(g,[-1 1]);
fplot(@(x) x,[-1 1]);
plot(xp,yp,'k');
hold off

