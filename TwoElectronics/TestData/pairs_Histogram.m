
A=load('init_undefined01311914.dat');
subplot(2,2,1);
histogram(A(:,3));
% hold on 
subplot(2,2,2);
histogram(A(:,9));
% %hold on
subplot(2,2,3);
 AA = load('z1z2.DAT');
 histogram(AA(:,1));

% hold on 

subplot(2,2,4);
histogram(AA(:,2));