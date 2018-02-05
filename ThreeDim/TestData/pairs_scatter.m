A=load('step_two_10w02060015.dat');

x1 = A(:,1);
y1 = A(:,2);
z1 = A(:,3);
px1 = A(:,4);
py1 = A(:,5);
pz1 = A(:,6);
x2 = A(:,7);
y2 = A(:,8);
z2 = A(:,9);
px2 = A(:,10);
py2 = A(:,11);
pz2 = A(:,12);


colorbar;
%  scatter(x1,x2,'.');
scatplot(pz1,pz2);