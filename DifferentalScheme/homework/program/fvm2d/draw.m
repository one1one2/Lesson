load('real40.dat');
load('weno40.dat');
load('real80.dat');
load('weno80.dat');
load('real160.dat');
load('weno160.dat');
load('real320.dat');
load('weno320.dat');

err40 = norm(real40-weno40,'fro')/40 
err80=norm(real80-weno80,'fro')/80 
err160=norm(real160-weno160,'fro')/160 
err320=norm(real320-weno320,'fro')/320 

log(err80/err40)/log(2)
log(err160/err80)/log(2)
log(err320/err160)/log(2)



