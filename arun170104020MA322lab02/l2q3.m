f = @(x)  x.^2 - x - 1;
ga = @(x) x.^2 - 1;
gb = @(x) 1 + 2*x - x.^2;
gc = @(x) (1 + 3*x - x.^2)/2;

% Part (a) With initial x =1
% x(1)=1;
% i=1;
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% while(abs(f(x))> 1e-2)
%     i = i+1;
%     x(i) = ga(x(i-1));
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end


% Part (a) With initial x =2
% x(1)=2;
% i=1;
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% while(abs(f(x)) > 1e-2)
%     i = i+1;
%     x(i) = ga(x(i-1));
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end

% Part (b) with initial x = 1
% x(1)=1;
% i=1;
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% while(abs(f(x))> 1e-2)
%     i = i+1;
%     x(i) = gb(x(i-1));
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end

% Part (b) with initial x = 2
% x(1)=2;
% i=1;
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% while(abs(f(x))> 1e-2)
%     i = i+1;
%     x(i) = gb(x(i-1));
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end


%Part (c) with initial x = 1

x(1)=1;
i=1;
fprintf('\n\t\t TABLE \n');
fprintf('\t n \t\t x_n \t\t f(x_n)\n');
fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
xp=[x(1)];
yp=[x(1)];
while(abs(f(x))> 1e-2)
    i = i+1;
    x(i) = gc(x(i-1));
    xp=[xp x(i-1)];
    xp=[xp x(i)];
    yp=[yp x(i)];
    yp=[yp x(i)];
    fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
end

hold on;
fplot(gc,[1 2]);
fplot(@(x) x,[1 2]);
plot(xp,yp,'k');
hold off

%Part (c) with initial x = 2
% x(1)=2;
% i=1;
% fprintf('\n\t\t TABLE \n');
% fprintf('\t n \t\t x_n \t\t f(x_n)\n');
% fprintf('\t %d \t\t %f \t %f\n',1,x(1),f(x(1)));
% xp=[x(1)];
% yp=[x(1)];
% while(abs(f(x))> 1e-2)
%     i = i+1;
%     x(i) = gc(x(i-1));
%     xp=[xp x(i-1)];
%     xp=[xp x(i)];
%     yp=[yp x(i)];
%     yp=[yp x(i)];
%     fprintf('\t %d \t\t %f \t %f\n',i,x(i),f(x(i)));
% end
% 
% hold on;
% fplot(gc,[1 2]);
% fplot(@(x) x,[1 2]);
% plot(xp,yp,'k');
% hold off