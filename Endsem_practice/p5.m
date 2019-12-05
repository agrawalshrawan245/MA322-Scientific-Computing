clear all;
clc;

% f = @(t,y,lambda) (lambda*y + (1-lambda)*cos(t) - (1+lambda)*sin(t));
% exact_sol = @(t) (sin(t) + cos(t));
% 
% ti = 0;
% tf = 10;
% h = 0.01;
% t = ti:h:tf;
% n = length(t);
% y(1) = 1; 
%Explicit Euler yn+1 =  yn + h*fn, using lambda = -5
% for i = 2:n
%     y(i) = y(i-1) + h*f(t(i-1), y(i-1), -5);
% end
% 
% plot(t,exact_sol(t));
% hold on;
% plot(t,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');

%Implicit Euler using yn+1 = yn + h*f(tn+1, yn+1)
% for i = 2:n
%     y_old = y(i-1) + h*f(t(i), y(i-1) + h*f(t(i-1), y(i-1), 5), 5);
%     y_new = y(i-1) + h*f(t(i), y_old, 5);
%     while abs(y_new - y_old) > 1e-5
%         save = y_new;
%         y_new = y(i-1) + h*f(t(i), y_new, 5);
%         y_old = save;
%     end
%     y(i) = y_new;
% end
% plot(t,exact_sol(t));
% hold on;
% plot(t,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');


% f = @(x,y) (x + y);
% exact_sol = @(x) (2*exp(x) -x -1);
% a = 0;
% b = 1;
% h = 0.005;
% x = a:h:b;
% n = length(x);
% y(1) = 1;
% 
% %Euler's method
% for i = 2:n
%     y(i) = y(i-1) + h*f(x(i-1), y(i-1));
% end
% 
% figure(1);
% plot(x,exact_sol(x));
% hold on;
% plot(x,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');
% 
% %Modified Euler's method
% z(1) = 1;
% for i = 2:n
%     z(i) = z(i-1) + (h/2)*(f(x(i-1), z(i-1)) + f(x(i), z(i-1) + h*f(x(i-1), z(i-1))));
% end
% 
% figure(2);
% plot(x,exact_sol(x));
% hold on;
% plot(x,z,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');


% f = @(t,y,lambda) (lambda*y + (1/(1 + t^2)) - lambda*atan(t));
% exact_sol = @(t) atan(t);
% 
% ti = 0;
% tf = 1;
% h = 0.001;
% t = ti:h:tf;
% n = length(t);
% y(1) = 0;
% lambda = -1;
% %Explicit Euler's method
% for i = 2:n
%     y(i) = y(i-1) + h*f(t(i-1), y(i-1), lambda);
% end
% 
% figure(1);
% plot(t,exact_sol(t));
% hold on;
% plot(t,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');
% 
% %Implcit Euler's method
% for i = 2:n
%     y_old = y(i-1) + h*f(t(i), y(i-1) + h*f(t(i-1), y(i-1), lambda), lambda);
%     y_new = y(i-1) + h*(f(t(i), y_old, lambda));
%     while abs(y_old - y_new) > 1e-5
%         y_old = y_new;
%         y_new = y(i-1) + h*f(t(i), y_old, lambda); 
%     end
%     y(i) = y_new;
% end
% 
% figure(2);
% plot(t,exact_sol(t));
% hold on;
% plot(t,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');
% 
% %Trapezoidal method
% for i = 2:n
%     y_old = y(i-1) + (h/2)*(f(t(i-1), y(i-1), lambda) + f(t(i), y(i-1) + h*f(t(i-1), y(i-1), lambda), lambda));
%     y_new = y(i-1) + (h/2)*(f(t(i-1), y(i-1), lambda) + f(t(i), y_old, lambda));
%     while abs(y_old - y_new) > 1e-5
%         y_old  = y_new;
%         y_new = y(i-1) + (h/2)*(f(t(i-1), y(i-1), lambda) + f(t(i), y_old, lambda)); 
%     end
%     y(i) = y_new;
% end
% 
% figure(3);
% plot(t,exact_sol(t));
% hold on;
% plot(t,y,'--', 'LineWidth', 2);
% legend('Exact solution', 'Approximation');

f = @(x) (log(x)*x^2);
a = 1;
b = 1.5;
n = 5;
changed_f = change_range(f, a, b);
approx = Gaussian_quad(changed_f, a, b, n);
function approx = Gaussian_quad(f, a, b, n)
    approx = 0;
    syms x;
    roots = vpasolve(legendreP(n,x));
    roots = double(reshape(roots,1,n));
    der = polyder(sym2poly(legendreP(n,x)));
    for i = 1:n
        w(i) = 2/((1 - roots(i)^2)*(polyval(der, roots(i)))^2);
        approx = approx + w(i)*f(roots(i));
    end
end

function f = change_range(f, a, b)
    syms t;
    outer = (b-a)/2;
    inner = ((b-a)*t + (b+a))/2;
    g = outer*f(inner);
    f = matlabFunction(g);
end