load('weno40.dat');
load('weno80.dat');
load('weno160.dat');
load('weno320.dat');

N = 320;
dat = weno320;

dat40 = dat;
for i=1:N
  for j=1:N
    dat40(i,j) = weno40(fix((i-1)/(N/40))+1,fix((j-1)/(N/40)+1));
  end
end 


dat80 = dat;
for i=1:N
  for j=1:N
    dat80(i,j) = weno80(fix((i-1)/(N/80))+1,fix((j-1)/(N/80)+1));
  end
end


dat160 = dat;
for i=1:N
  for j=1:N
    dat160(i,j) = weno160(fix((i-1)/(N/160))+1,fix((j-1)/(N/160)+1));
  end
end 


err1=norm(dat40-dat,'fro')/320
err2=norm(dat80-dat,'fro')/320
err3=norm(dat160-dat,'fro')/320

log(err2/err1)/log(2)
log(err3/err2)/log(2)
