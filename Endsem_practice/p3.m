clear all;
clc;

% f = @(x) (x^3 - x -1);
% g = @(x) (x + 1).^(1/3);
% tol = 1e-2;
% x(1) = 1;
% i = 1;
% 
% xx = 1:0.01:2;
% plot(xx,g(xx));
% hold on;
% plot(xx,xx);
% xl = [x(1)];
% yl = [g(x(1))];
% while abs(f(x(i))) > tol
%     i = i + 1;
%     x(i) = g(x(i-1));
%     xl = [xl x(i) x(i)];
%     yl = [yl x(i) g(x(i))];
% end
% 
% hold on;
% plot(xl, yl);


% f = @(x) (3.*x.*x - exp(x));
% g = @(x) (exp(x)/3).^(1/2);
% 
% xx = 0.8:0.01:1;
% plot(xx, g(xx));
% hold on;
% plot(xx, xx);
% 
% x(1) = 1;
% xl = [x(1)];
% yl = [g(x(1))];
% 
% i = 1;
% while abs(f(x(i))) > 1e-3
%     i = i + 1;
%     x(i) = g(x(i-1));
%     xl = [xl x(i) x(i)];
%     yl = [yl x(i) g(x(i))]; 
% end
% 
% hold on;
% plot(xl ,yl ,'m');


% f = @(x) (x - cos(x));
% g = @(x) (cos(x));
% 
% x(1) = 0.9;
% xl = [x(1)];
% yl = [g(x(1))];
% i = 1;
% 
% xx = 0.5:0.01:1;
% plot(xx, g(xx));
% hold on;
% plot(xx,xx);
% 
% while abs(f(x)) > 1e-4
%     i = i + 1;
%     x(i) = g(x(i-1));
%     xl = [xl x(i) x(i)];
%     yl = [yl x(i) g(x(i))];
% end
% 
% hold on;
% plot(xl,yl,'m');


% f = @(x) (x.^2 - 0.75);
% g = @(x) (x.^2 + x -0.75);
% x(1) = 0.8;
% xl = [x(1)];
% yl = [g(x(1))];
% i = 1;
% 
% xx = -1:0.01:1.7;
% plot(xx, g(xx));
% hold on;
% plot(xx,xx);
% 
% while abs(f(x)) > 1e-4
%     i =  i + 1;
%     x(i) = g(x(i-1));
%     xl = [xl x(i) x(i)];
%     yl = [yl x(i) g(x(i))];
% end
% 
% hold on;
% plot(xl, yl, 'm');


% x(1) = 1;
% x(2) = 1.2;
% x(3) = 1.4;
% 
% f = @(x) (x.^3 - x - 2);
% fdiff = @(x,y) ((f(x) - f(y))/ (x-y));
% fdiff2 = @(x,y,z) ((fdiff(x,y) - fdiff(y,z))/(x-z));
% 
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% fprintf('\t %d \t\t %f \t %f\n',2,x(2),f(x(2)));
% fprintf('\t %d \t\t %f \t %f\n',3,x(2),f(x(3)));
% i = 3;
% while(abs(f(x(i))) > 1e-3 )
%     i = i + 1;
%     w = fdiff(x(i-1), x(i-2)) + fdiff(x(i-1), x(i-3)) - fdiff(x(i-2), x(i-3));
%     x(i) = x(i-1)-(2*f(x(i-1))/(max(w+sqrt(w^2-4*f(x(i-1))*fdiff2(x(i-1),x(i-2),x(i-3))),w-sqrt(w^2-4*f(x(i-1))*fdiff2(x(i-1),x(i-2),x(i-3))))));       
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end
% 
% root = x(i)  


% y(1) = 1.5;
% y(2) = 1.4;
% y(3) = 1.3;
% 
% f = @(x) (1 + 2*x - tan(x));
% fdiff = @(x,y) ((f(x) - f(y))/ (x-y));
% fdiff2 = @(x,y,z) ((fdiff(x,y) - fdiff(y,z))/(x-z));
% 
% i=3;
% 
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(1)));
% fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(2)));
% fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(3)));
% 
% while(abs(f(y(i))) > 1e-3)
%     i = i+1;
%     w = fdiff(y(i-1), y(i-2)) + fdiff(y(i-1), y(i-3)) - fdiff(y(i-2), y(i-3));
%     y(i) = y(i-1)-(2*f(y(i-1))/(max(w+sqrt(w^2-4*f(y(i-1))*fdiff2(y(i-1),y(i-2),y(i-3))),w-sqrt(w^2-4*f(y(i-1))*fdiff2(y(i-1),y(i-2),y(i-3))))));
%     fprintf('\t %d \t\t %f \t %f\n',i,y(i),f(y(i)));
% end
% 
% root = y(i)

f = @(x) (x^3 - x - 2);
fdiff1 = @(x,y) (f(x) - f(y))/(x-y);
fdiff2 = @(x,y,z) (fdiff1(x,y) - fdiff1(y,z))/(x-z);

x(1) = 1;
x(2) = 1.2;
x(3) = 1.4;
i = 3;
while abs(f(x(i))) > 1e-4
    i = i + 1;
    w = fdiff1(x(i-1), x(i-2))  + fdiff1(x(i-1), x(i-3)) - fdiff1(x(i-2), x(i-3));
    x(i) = x(i-1) - (2*f(x(i-1)))/max((w + sqrt(w^2 - 4*f(x(i-1))*fdiff2(x(i-1), x(i-2), x(i-3)))),(w - sqrt(w^2 - 4*f(x(i-1))*fdiff2(x(i-1), x(i-2), x(i-3)))));
end

root =  x(i)