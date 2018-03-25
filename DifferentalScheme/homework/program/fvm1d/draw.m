
a_LF=load('LF.dat');
a_LW=load('LW.dat');
a_Go=load('Go.dat');
x=linspace(-1,1,size(a_LF,2));
plot(x,a_LF,'o',x,a_LW,'o',x,a_Go,'o');
legend('LF','LW','Godunov');
