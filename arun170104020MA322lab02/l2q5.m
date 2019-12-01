x(1) = 1;
x(2) = 1.2;
x(3) = 1.4;

f = @(x) (x.^3 - x - 2);
fdiff = @(x,y) ((f(x) - f(y))/ (x-y));
fdiff2 = @(x,y,z) ((fdiff(x,y) - fdiff(y,z))/(x-z));

fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(2)));
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(3)));
i = 3;
while(abs(f(x(i))) > 1e-3 )
    i = i + 1;
    w = fdiff(x(i-1), x(i-2)) + fdiff(x(i-1), x(i-3)) - fdiff(x(i-2), x(i-3));
    x(i) = x(i-1)-(2*f(x(i-1))/(max(w+sqrt(w^2-4*f(x(i-1))*fdiff2(x(i-1),x(i-2),x(i-3))),w-sqrt(w^2-4*f(x(i-1))*fdiff2(x(i-1),x(i-2),x(i-3))))));       
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
end

root = x(i)  


y(1) = 1.5;
y(2) = 1.4;
y(3) = 1.3;

f = @(x) (1 + 2*x - tan(x));
fdiff = @(x,y) ((f(x) - f(y))/ (x-y));
fdiff2 = @(x,y,z) ((fdiff(x,y) - fdiff(y,z))/(x-z));

i=3;

fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(1)));
fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(2)));
fprintf('\t %d \t\t %f \t %f\n',1,y(1),f(y(3)));

while(abs(f(y(i))) > 1e-3)
    i = i+1;
    w = fdiff(y(i-1), y(i-2)) + fdiff(y(i-1), y(i-3)) - fdiff(y(i-2), y(i-3));
    y(i) = y(i-1)-(2*f(y(i-1))/(max(w+sqrt(w^2-4*f(y(i-1))*fdiff2(y(i-1),y(i-2),y(i-3))),w-sqrt(w^2-4*f(y(i-1))*fdiff2(y(i-1),y(i-2),y(i-3))))));
    fprintf('\t %d \t\t %f \t %f\n',i,y(i),f(y(i)));
end

root = y(i)


