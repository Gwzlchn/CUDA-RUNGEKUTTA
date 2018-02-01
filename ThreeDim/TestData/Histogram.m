A=load('init_1w.dat');

histogram(A(:,3))
hold on
histogram(A(:,9))

%axis([-5,5,-5,5]);