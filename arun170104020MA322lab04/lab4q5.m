% Part a
x = [0.1 0.2 0.3 0.4];
fx = [-0.29004986 -0.56079734 -0.81401972 -1.0526302];

ans =0;

for i=1:4
    tmp = 1;
    for j =1:4
        if i~=j       
            tmp = tmp*( (0.18-x(j))/(x(i)-x(j)) );
        end
    end
    ans = ans + tmp*fx(i);
end

parta = ans

%Part b

x = [-1 -0.5 0 0.5];
gx = [0.86199480 0.95802009 1.0986123 1.2943767];

ans =0;

for i=1:4
    tmp = 1;
    for j =1:4
        if i~=j       
            tmp = tmp*( (0.25-x(j))/(x(i)-x(j)) );
        end
    end
    ans = ans + tmp*gx(i);
end

partb = ans


