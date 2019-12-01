f = @(x) x^5 -2*x^4 -2*x^3 + 8*x^2  - 7*x +2;
syms x;
F = x^5 -2*x^4 -2*x^3 + 8*x^2  - 7*x +2;
derf = matlabFunction(diff(F));
df = derf;

x = 1.3;
it = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');

while(abs(f(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
    it = it+1;
    x = x - 4*(f(x)/df(x));
end 
fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
modified = x
it


g = @(x) x^5 -8*x^4 +25*x^3 -38*x^2 +28*x -8;
syms x;
G = x^5 -8*x^4 +25*x^3 -38*x^2 +28*x -8;
dg = matlabFunction(diff(G));


x = 0;
it = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');

while(abs(g(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
    it = it+1;
    x = x - 2*(g(x)/dg(x));
end 
fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));

x = 3;
it = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');

while(abs(g(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
    it = it+1;
    x = x - 3*(g(x)/dg(x));
end 
fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
