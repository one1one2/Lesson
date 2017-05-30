a = load('1.dat');
rho = a(:,1);
u = a(:,2)./a(:,1);
gamma = 1.4;
p = (gamma - 1)*(a(:,3) - 0.5*rho.*u.*u);
x = linspace(-5,5,size(rho,1));

plot(x,rho)
