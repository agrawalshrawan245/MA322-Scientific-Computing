%To run uncomment the respective section and do clear all and run to see
%the output

%%% PROBLEM 1

% lower = -4;
% upper = -3;
% mid(1) = (lower + upper)/2;
% i=1;
% while abs(func(mid(i))) > 0.001
%    i=i+1;
%    mid(i)= (lower + upper)/2;
%    if (func(lower) * func(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% for i=3:9
%     e(i-2) = mid(i) - mid(i-1);
% end
% semilogx([1:7], e);
% 
% for i=1:7
%     er(i) = log2(abs(e(i+1)/e(i)));
% end
% 
% function f = func(n)
%     f = exp(n)-sin(n);
% end



%%% PROBLEM 2

% % For part 1
% lower = 0.5;
% upper = 1.5;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f1(mid(i))) > 0.001
%    i= i+1;
%    mid(i) = (lower + upper)/2;
%    if (f1(lower) * f1(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 1 is x = ' + string(mid(i)));
% 
% function f = f1(x)
%     f = 2 + cos(exp(x)-2)-exp(x);
% end


% % For Part 2

% lower = 0;
% upper = 1;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f2(mid(i))) > 0.001
%     i = i+1;
%    mid(i) = (lower + upper)/2;
%    if (f2(lower) * f2(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 2 is x = ' +string(mid(i)));
% 
% function f = f2(x)
%     f = x - power(2, -x);
% end


% % For part 3

% lower = 0;
% upper = 1;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f3(mid(i))) > 0.001
%     i = i+1;
%    mid(i) = (lower + upper)/2;
%    if (f3(lower) * f3(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 3 is x = ' +string(mid(i)));
% 
% function f = f3(x)
%     f = exp(x) - power(x,2) + 3*x -2; 
% end

% % For Part 4

% -3<=x<=-2
% lower = -3;
% upper = -2;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f4(mid(i))) > 0.001
%    mid(i) = (lower + upper)/2;
%    if (f4(lower) * f4(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end

% disp('Solution for function 4 is x for range -3 to -2 is = ' +string(mid(i)));

%  -1<=x<=0
% lower = -1;
% upper = -0;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f4(mid(i))) > 0.001
%    mid(i) = (lower + upper)/2;
%    if (f4(lower) * f4(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 4 for x in range -1 to 0 is = ' +string(mid(i)));
% 
% function f = f4(x)
%     f = 2*x*cos(2*x) - (x+1)*(x+1);
% end


% % For Part 5

% % 0.2<=x<=0.3
% lower = 0.2;
% upper = 0.3;
% i=1;
% mid(i) = (lower + upper)/2;
% while abs(f5(mid(i))) > 0.001
%    mid(i) = (lower + upper)/2;
%    if (f5(lower) * f5(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 5 for x in range 0.2 to 0.3 is = ' +string(mid(i)));
% 
% % 1.2<=x<=1.3
% lower = 1.2;
% upper = 1.3;
% i=1;
% mid(i)= (lower + upper)/2;
% while abs(f5(mid(i))) > 0.001
%    mid = (lower + upper)/2;
%    if (f5(lower) * f5(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% 
% disp('Solution for function 5 for x in range 1.2 to 1.3 is = ' +string(mid(i)));
% 
% function f = f5(x)
%     f = x*cos(x) - 2*x*x + 3*x -1;
% end


% % PROBLEM  3

% lower = 2.5;
% upper = 3;
% i=1;
% mid(i) = (lower + upper)/2;
% 
% while abs(func(mid(i))) > 0.0001
%     i = i+1;
%    mid(i) = (lower + upper)/2;
%    val(i-1) = func(mid(i));
%    if (func(lower) * func(mid(i))) < 0
%        upper = mid(i);
%    else 
%        lower=mid(i);
%    end
% end
% val = val.';
% mid = mid.';
% function f = func(x)
%     f = power(x, 3) - 25;
% end
% 

% % PROBLEM  4

% x=1;
% i=0;
% while(abs(f(x))>1e-10)
%     i = i+1;
%     if i==5
%         disp('Value of f(x) for fourth iterate for part (i) = ' + string(f(x)));
%     end
%     x = x*power((1 + (7 - power(x,5))/(x*x)),3);
% end
% disp('Number of iterations in part (i) = ' + string(i));
% 
% x=1;
% i=0;
% while(abs(f(x))>1e-10)
%     i = i+1;
%     if i==5
%         disp('Value of f(x) for fourth iterate for part(ii) = ' + string(f(x)));
%     end
%     x = x - ((power(x,5)-7)/(x*x));
% end
% disp('Number of iterations in part (ii) = ' + string(i));
% 
% x=1;
% i=0;
% while(abs(f(x))>1e-10)
%     i = i+1;
%     if i==5
%         disp('Value of f(x) for fourth iterate for part(iii) = ' + string(f(x)));
%     end
%     x = x - ((power(x,5)-7)/(5*power(x,4)));
% end
% 
% disp('Number of iterations in part (iii) = ' + string(i) + ' and approximation is ' + string(x));
% 
% x=1;
% i=0;
% while(abs(f(x))>1e-10)
%     i = i+1;
%     if i==5
%         disp('Value of f(x) for fourth iterate for part(iv) = ' + string(f(x)));
%     end
%     x = x - ((power(x,5)-7)/(12));
% end
% disp('Number of iterations in part (iv) = ' + string(i) + ' and approximation is ' + string(x));
% 
% function out = f(x)
% out = power(x,5)-7;
% end


% % PROBLEM 5

% x = 1.7;
% while(abs(func1(x))>0.00001)
%     x = x - func1(x)/derfunc1(x);
% end
% 
% disp('Solution for function 1 is ' + string(x));
% 
% 
% x = 1.5;
% while(abs(func2(x))>0.00001)
%     x = x - func2(x)/derfunc2(x);
% end
% 
% disp('Solution for function 2 is ' + string(x));
% 
% 
% 
% x = 2.5;
% while(abs(func3(x))>0.00001)
%     x = x - func3(x)/derfunc3(x);
% end
% 
% disp('Solution for function 3 in range 2 to 3 is ' + string(x));
% 
% 
% x = 3.5;
% while(abs(func3(x))>0.00001)
%     x = x - func3(x)/derfunc3(x);
% end
% 
% disp('Solution for function 3 in range 3 to 4 is ' + string(x));
% 
% 
% x = 1.5;
% while(abs(func4(x))>0.00001)
%     x = x - func4(x)/derfunc4(x);
% end
% 
% disp('Solution for function 4 in range 1 to 2 is ' + string(x));
% 
% x = (exp(1) + 4)/2;
% while(abs(func4(x))>0.00001)
%     x = x - func4(x)/derfunc4(x);
% end
% 
% disp('Solution for function 4 in range e to 4 is ' + string(x));
% 
% x = 0.5;
% while(abs(func5(x))>0.00001)
%     x = x - func5(x)/derfunc5(x);
% end
% 
% disp('Solution for function 5 in range 0 to 1 is ' + string(x));
% 
% x = 4;
% while(abs(func5(x))>0.00001)
%     x = x - func5(x)/derfunc5(x);
% end
% 
% disp('Solution for function 5 in range 3 to 5 is ' + string(x));
% 
% x = 0.5;
% while(abs(func6(x))>0.00001)
%     x = x - func6(x)/derfunc6(x);
% end
% 
% disp('Solution for function 6 in range 0 to 1 is ' + string(x));
% 
% 
% x = 3.5;
% while(abs(func6(x))>0.00001)
%     x = x - func6(x)/derfunc6(x);
% end
% 
% disp('Solution for function 6 in range 3 to 4 is ' + string(x));
% 
% x = 6.5;
% while(abs(func6(x))>0.00001)
%     x = x - func6(x)/derfunc6(x);
% end
% 
% disp('Solution for function 6 in range 6 to 7 is ' + string(x));
% 
% function f = func1(x)
%     f = exp(x) + power(2,-x) + 2*cos(x) -6; 
% end
% 
% function f = derfunc1(x)
%     f = exp(x) - power(2,-x)*log(2) - 2*sin(x);
% end
% 
% function f = func2(x)
%     f = log(x-1) + cos(x-1); 
% end
% 
% function f = derfunc2(x)
%     f = (1/(x-1)) - sin(x-1);
% end
% 
% 
% function f = func3(x)
%     f = 2*x*cos(2*x) - power(x-2,2); 
% 
% end
% 
% function f = derfunc3(x)
%     f = 2*cos(2*x) - 4*x*sin(2*x) - 2*x +4;
% end
% 
% 
% function f = func4(x)
%     f = power(x-2,2)  - log(x); 
% end
% 
% function f = derfunc4(x)
%     f = 2*(x-2) - (1/x);
% end
% 
% function f = func5(x)
%     f = exp(x) - 3*x*x; 
% end
% 
% function f = derfunc5(x)
%     f = exp(x) - 6*x;
% end
% 
% function f = func6(x)
%     f = sin(x) - (1/exp(x)); 
% end
% 
% function f = derfunc6(x)
%     f = cos(x) + (1/exp(x));
% end


% % PROBLEM 6 IS IN SEPARATE FILE.


% % PROBLEM 7 has two files named Newton7 and Secant7 which can take
% required custom inputs.

% % PROBLEM 8 also has separate code files





