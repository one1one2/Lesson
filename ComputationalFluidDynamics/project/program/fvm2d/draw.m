rho=load('rho.dat');
x=linspace(0,4,size(rho,2));
y=linspace(0,1,size(rho,1));
contour(x,y,rho, 50)
axis equal
axis([0,3,0,1])
