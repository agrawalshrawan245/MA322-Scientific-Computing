clear all;
clc;

% p = @(x) -2;
% q = @(x) (-1);
% r = @(x) (-x);
% 
% a = 0;
% b = 1;
% h = 0.01;
% x = a:h:b;
% n = length(x);
% y(1) = 0;
% y(n) = 0;

%Forward difference
% A(1,1) = (2/(h^2)) - p(x(2))/h + q(x(2));
% A(1,2) = (-1/h^2) + p(x(2))/h;
% b(1) = r(x(2)) + (1/(h^2))*y(1);
% for i = 2:n-3
%     A(i,i-1) = (-1/(h^2));
%     A(i,i) = (2/(h^2)) - p(x(i+1))/h + q(x(i+1));
%     A(i,i+1) = (-1/(h^2)) + p(x(i+1))/h;
%     b(i) = r(x(i+1));
% end
% A(n-2, n-3) = (-1/(h^2));
% A(n-2, n-2) = (2/(h^2)) - p(x(n-1))/h + q(x(n-1));
% b(n-2) = r(x(n-1)) - ((-1/(h^2)) + p(x(n-1))/h)*y(n); 
% 
% u = A\(b');
% solution = [y(1) u' y(n)];
% plot(x,solution);

%Backward difference
% A(1,1) = (2/(h^2)) + p(x(2))/h + q(x(2));
% A(1,2) = (-1/h^2);
% b(1) = r(x(2)) + ((1/(h^2)) - p(x(2))/h)*y(1);
% for i = 2:n-3
%     A(i,i-1) = (-1/(h^2)) - p(x(i+1))/h;
%     A(i,i) = (2/(h^2)) + p(x(i+1))/h + q(x(i+1));
%     A(i,i+1) = (-1/(h^2));
%     b(i) = r(x(i+1));
% end
% A(n-2, n-3) = (-1/(h^2)) - p(x(n-1))/h;
% A(n-2, n-2) = (2/(h^2)) + p(x(n-1))/h + q(x(n-1));
% b(n-2) = r(x(n-1)) + (1/(h^2))*y(n); 
% 
% u = A\(b');
% solution = [y(1) u' y(n)];
% plot(x,solution);

%Central difference
% A(1,1) = (2/(h^2)) + q(x(2));
% A(1,2) = (-1/h^2) + p(x(2))/(2*h);
% b(1) = r(x(2)) + ((1/(h^2)) + p(x(2))/(2*h))*y(1);
% for i = 2:n-3
%     A(i,i-1) = (-1/(h^2)) - (p(x(i+1))/(2*h));
%     A(i,i) = (2/(h^2)) + q(x(i+1));
%     A(i,i+1) = (-1/(h^2)) + p(x(i+1))/(2*h);
%     b(i) = r(x(i+1));
% end
% A(n-2, n-3) = (-1/(h^2))- (p(x(n-1))/(2*h));
% A(n-2, n-2) = (2/(h^2)) + q(x(n-1));
% b(n-2) = r(x(n-1)) - ((-1/(h^2)) + p(x(n-1))/(2*h))*y(n); 
% 
% u = A\(b');
% solution = [y(1) u' y(n)];
% plot(x,solution);

% p = @(x) 3;
% q = @(x) -2;
% r = @(x) -(2);
% a = 0;
% b = 1;
% h = 0.01;
% x = a:h:b;
% n = length(x);
% 
% b(1) = h;
% A(1,1) = -1-h;
% A(1,2) = 1;
% for i = 2:n-1
%     A(i,i-1) = (-1/(h^2)) - (p(x(i))/(2*h));
%     A(i,i) = (2/(h^2)) + q(x(i));
%     A(i,i+1) = (-1/(h^2)) + p(x(i))/(2*h);
%     b(i) = r(x(i));
% end
% b(n) = h;
% A(n, n-1) = -1;
% A(n, n) = 1 + h;
% 
% u = A\(b');
% plot(x,u);
% hold on;
% plot(x,exact_sol(x));



