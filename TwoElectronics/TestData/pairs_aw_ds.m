A=load('DS_10w02060041.dat');

%scatter(A(:,3),A(:,9))
AA = load('AW_10w02060041.dat');

PI = 3.1415926;
omega = 0.057;
t0 = 2 * PI / omega;
h=0.027;
B =zeros(80000,1);
for i = 1:80000
    B(i) = i*h*0.5/t0;
end
% scatter(B,A,'.')
% hold on 
% scatter(B,AA,'.')
% save('time.txt','B')

 AAA=load('1E(t).dat');
 scatter(AAA(:,1),AAA(:,2),'.');
 hold on;
 scatter(AAA(:,1),AAA(:,3),'.');

scatter(AAA(:,1),-A,'.')
hold on 
scatter(AAA(:,1),AA,'.')
