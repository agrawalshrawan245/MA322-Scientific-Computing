clear all;
clc;

fa = @(x) (x.^2).*log(x);
fb = @(x) (exp(3*x).*sin(2*x));
fc = @(x) (2/(x.^2 - 4));
fd = @(x) (2.*x/(x.^2 - 4));

%Part a
a = 1;
b = 1.5;
fa = change_range(fa,a,b);
approx = Gaussian_quadrature(fa,a,b,2);
fprintf('Approximation with Gaussian quadrature using n=2 is %.8f.\n',approx);
approx = Gaussian_quadrature(fa,a,b,3);
fprintf('Approximation with Gaussian quadrature using n=3 is %.8f.\n',approx);
approx = Gaussian_quadrature(fa,a,b,4);
fprintf('Approximation with Gaussian quadrature using n=4 is %.8f.\n',approx);
approx = Gaussian_quadrature(fa,a,b,5);
fprintf('Approximation with Gaussian quadrature using n=5 is %.8f.\n\n',approx);

%Part b
a = 0;
b = pi/4;
fb = change_range(fb,a,b);
approx = Gaussian_quadrature(fb,a,b,2);
fprintf('Approximation with Gaussian quadrature using n=2 is %.8f.\n',approx);
approx = Gaussian_quadrature(fb,a,b,3);
fprintf('Approximation with Gaussian quadrature using n=3 is %.8f.\n',approx);
approx = Gaussian_quadrature(fb,a,b,4);
fprintf('Approximation with Gaussian quadrature using n=4 is %.8f.\n',approx);
approx = Gaussian_quadrature(fb,a,b,5);
fprintf('Approximation with Gaussian quadrature using n=5 is %.8f.\n\n',approx);

%Part c
a = 0;
b = 0.35;
fc = change_range(fc,a,b);
approx = Gaussian_quadrature(fc,a,b,2);
fprintf('Approximation with Gaussian quadrature using n=2 is %.8f.\n',approx);
approx = Gaussian_quadrature(fc,a,b,3);
fprintf('Approximation with Gaussian quadrature using n=3 is %.8f.\n',approx);
approx = Gaussian_quadrature(fc,a,b,4);
fprintf('Approximation with Gaussian quadrature using n=4 is %.8f.\n',approx);
approx = Gaussian_quadrature(fc,a,b,5);
fprintf('Approximation with Gaussian quadrature using n=5 is %.8f.\n\n',approx);

%Part d
a = 1;
b = 1.6;
fd = change_range(fd,a,b);
approx = Gaussian_quadrature(fd,a,b,2);
fprintf('Approximation with Gaussian quadrature using n=2 is %.8f.\n',approx);
approx = Gaussian_quadrature(fd,a,b,3);
fprintf('Approximation with Gaussian quadrature using n=3 is %.8f.\n',approx);
approx = Gaussian_quadrature(fd,a,b,4);
fprintf('Approximation with Gaussian quadrature using n=4 is %.8f.\n',approx);
approx = Gaussian_quadrature(fd,a,b,5);
fprintf('Approximation with Gaussian quadrature using n=5 is %.8f.\n',approx);

function approx = Gaussian_quadrature(f,a,b,n)
    approx = 0;
    syms x
    roots = vpasolve(legendreP(n,x));
    roots = double(reshape(roots,1,n));
    der = polyder(sym2poly(legendreP(n,x)));
    w = [];
    for i = 1:n
        w(i) = 2./((1 - roots(i).^2).*((polyval(der,roots(i))).^2));
        approx = approx + w(i)*f(roots(i));
    end
end

function f = change_range(F,a,b)
   syms t;
   outer = (b-a)/2;
   inner = ((b-a)*t + (b+a))/2;
   f = outer*F(inner);
   f = matlabFunction(f);
end