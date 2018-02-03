A=load('DS_100002040328.dat');

%scatter(A(:,3),A(:,9))
AA = load('AW_1000.dat');

PI = 3.1415926;
omega = 0.057;
t0 = 2 * PI / omega;
h=0.027;
B =zeros(80000,1);
for i = 1:80000
    B(i) = i*h*0.5/t0;
end
scatter(B,A)
hold on 
scatter(B,AA)
% save('time.txt','B')