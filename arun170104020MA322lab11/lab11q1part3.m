clear all;
clc;
clf;

h = 0.02;
xa = 0;
xb = 1;
ya = 0;
yb = 1;
x = [xa:h:xb];
y = [ya:h:yb];
nx = (xb-xa)/h+1; %no of nodal points in x-axis
ny = (yb-ya)/h+1;
n = nx*ny;


%Part c
exact_sol = @(x,y) (2*x - y);
f = @(x,y) (2*x - y);
bc_ya = @(x) (2*x);
bc_yb = @(x) (2*x - 1);
bc_xa = @(h, y) (h*(2-y)); 
bc_xb = @(y) (2 - y);



k = 1;
for i=2:nx-1
    for j=2:ny-1
    cur=i + (j-1)*nx;
    %A(cur,cur)=-4;
    first_index(k) = cur;
    second_index(k) = cur;
    s(k) = (h^2 - 4);
    k = k + 1;
    r(cur, 1) = (h^2)*f(x(i), y(j));
    %A(cur,cur-1)=1;
    first_index(k) = cur;
    second_index(k) = cur - 1;
    s(k) = 1;
    k = k + 1;
    %A(cur,cur+1)=1;
    first_index(k) = cur;
    second_index(k) = cur + 1;
    s(k) = 1;
    k = k + 1;
    %A(cur,cur-nx)=1;
    first_index(k) = cur;
    second_index(k) = cur - nx;
    s(k) = 1;
    k = k + 1;
    %A(cur,cur+nx)=1;
    first_index(k) = cur;
    second_index(k) = cur + nx;
    s(k) = 1;
    k = k + 1;
    end
end

%x=0 bc
for j=1:ny
  cur=1+(j-1)*nx;
  %A(cur,cur)=1;
  first_index(k) = cur;
  second_index(k) = cur;
  s(k) =h - 1;
  k = k + 1;
  first_index(k) = cur;
  second_index(k) = cur+1;
  s(k) = 1;
  k = k + 1;
  r(cur,1)=bc_xa(h,y(j));
end

%x=1 bc
for j=1:ny
  cur=nx+(j-1)*nx;
  %A(cur,cur)=1;
  first_index(k) = cur;
  second_index(k) = cur;
  s(k) = 1;
  k = k + 1;
  r(cur,1)=bc_xb(y(j));
end

%y=0 bc
for i=2:nx-1
  cur= i;
  %A(cur,cur)=1;
  first_index(k) = cur;
  second_index(k) = cur;
  s(k) = 1;
  k = k + 1;
  r(cur,1)=bc_ya(x(i));
end

%y=1 bc
for i=2:nx-1
  cur=i+(ny-1)*nx;
  %A(cur,cur)=1;
  first_index(k) = cur;
  second_index(k) = cur;
  s(k) = 1;
  k = k + 1;
  r(cur,1)=bc_yb(x(i));
end

S = sparse(first_index, second_index, s, n, n);
ans = S\r;


for i=1:nx
    for j=1:ny
    cur=i + (j-1)*nx;
    u(i,j) = ans(cur);
    exact(i,j) = exact_sol(x(i),y(j));
    err(i,j) = abs(u(i,j) - exact(i,j));
    end
end


figure(1);
[X,Y]=meshgrid(x,y);
surf(X,Y,u);
title('Surface plot of Numerical solution');
xlabel('x(i)');
ylabel('y(i)');
zlabel('u(x(i), y(i))');
figure(2);
surf(X,Y,exact);
title('Surface plot of Exact solution');
xlabel('x(i)');
ylabel('y(i)');
zlabel('u(x(i), y(i))');
figure(3);
surf(X,Y,err);
title('Surface plot of Error');
xlabel('x(i)');
ylabel('y(i)');
zlabel('e(x(i), y(i))');
figure(4);
contour(X,Y,u);
title('Contour plot of Numerical solution');
xlabel('x(i)');
ylabel('y(i)');
zlabel('u(x(i), y(i))');
figure(5);
contour(X,Y,exact);
title('Contour plot of Exact solution');
xlabel('x(i)');
ylabel('y(i)');
zlabel('u(x(i), y(i))');