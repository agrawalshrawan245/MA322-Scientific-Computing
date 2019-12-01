f = @(x) cos(x + sqrt(2)) + x*((x/2) + sqrt(2));
syms x;
F = cos(x + sqrt(2)) + x*((x/(2)) + sqrt(2));
derf = matlabFunction(diff(F));
df = derf;

x = -1;
itr = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
while(abs(f(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',itr,x,f(x));
    itr = itr + 1;
    x = x - (f(x)/df(x));
end 
fprintf('\t %d \t\t %f \t %f\n',itr,x,f(x));

p=1;
F = diff(F); 
while(abs(derf(x)) < 1e-3)
    p = p+1;
    F = diff(F); 
    derf = matlabFunction(F);
end

x = -1;
it = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
while(abs(f(x)) > 1e-5)
    fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));
    it = it+1;
    x = x - p*(f(x)/df(x));
end 
fprintf('\t %d \t\t %f \t %f\n',it,x,f(x));

g = @(x) exp(6*x) + 3*exp(2*x)*(log(2))^2 - log(8)*exp(4*x) - (log(2))^3;
syms x;
G= exp(6*x) + 3*exp(2*x)*(log(2))^2 - log(8)*exp(4*x) - (log(2))^3;
derg = matlabFunction(diff(G));
dg = derg;


x = -1;
itr = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
while(abs(g(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',itr,x,g(x));
    itr = itr +1;
    x = x - (g(x)/dg(x));
end 
fprintf('\t %d \t\t %f \t %f\n',itr,x,g(x));
p=1;
while(abs(derg(x)) < 1e-5)
    p = p+1;
    G = diff(G); 
    derg = matlabFunction(G);
end

x = -1;
it = 0;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
while(abs(g(x)) > 1e-15)
    fprintf('\t %d \t\t %f \t %f\n',it,x,g(x));
    it = it+1;
    x = x - p*(g(x)/dg(x));
end 
fprintf('\t %d \t\t %f \t %f\n',it,x,g(x));