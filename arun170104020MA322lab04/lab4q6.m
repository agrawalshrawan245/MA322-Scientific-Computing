clear all;
clc;
x = [1950 1960 1970 1980 1990 2000];
fx = [151326 179323 203302 226542 249633 281422];

size = 4;
mat = zeros(size,size+1);

for i = 1:size
    mat(i,1) = x(i);
    mat(i,2) = fx(i);
end

for j = 3:size+1
    for i = j-1:size
        mat(i,j) = (mat(i,j-1) - mat(i-1,j-1))/(mat(i,1) - mat(i-j+2,1));
    end
end

year = 2020;
pre = [(year-x(1))];
for i=2:size-1
    pre(i) = pre(i-1)*(year-x(i));
end

population = fx(1);
for i =1:size-1
    population = population + pre(i)*mat(i+1,i+2);
end

fprintf('population = %f\n', population);