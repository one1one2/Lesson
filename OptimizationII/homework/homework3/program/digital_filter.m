function [A F] = digital_filter(x)
phi = zeros(41,1);
phi(1:6) = 0:0.01:0.05;
phi(7:20) = 0.07:0.03:0.46;
phi(21:22) = 0.5:0.04:0.54;
phi(23:24) = 0.57:0.1:0.67;
phi(25:35) = 0.63:0.03:0.93;
phi(36:41) = 0.95:0.01:1;
for i=1:41
  A{i} = Fun(x, phi(i));
end 
F = A{1}.f;
for i = 2:41
  F = max(F, A{i}.f);
end 

end 


function C = Fun(x, phi)
theta = pi*phi;
n = length(x);
k = (n-1)/4;
a = x(1:4:n-4);
b = x(2:4:n-3);
c = x(3:4:n-2);
d = x(4:4:n-1);

H = x(n);
for i = 1:k
  H = H*(1+a(i)^2+b(i)^2+2*b(i)*(2*cos(theta)^2-1)+2*a(i)*(1+b(i))*cos(theta))^(0.5); 
  H = H/(1+c(i)^2+d(i)^2+2*d(i)*(2*cos(theta)^2-1)+2*c(i)*(1+d(i))*cos(theta))^(0.5); 
end 
if H >= abs(1-2*phi)
  flag = 1;
else 
  flag = -1;
end 

C.f = abs(H - abs(1-2*phi));

C.g = zeros(n,1);

for i = 1:k
  K1 = 1+a(i)^2+b(i)^2+2*b(i)*(2*cos(theta)^2-1)+2*a(i)*(1+b(i))*cos(theta); 
  K2 = 1+c(i)^2+d(i)^2+2*d(i)*(2*cos(theta)^2-1)+2*c(i)*(1+d(i))*cos(theta);
  C.g(4*i-3) = 0.5*H/K1*(2*a(i)+2*(1+b(i))*cos(theta));
  C.g(4*i-2) = 0.5*H/K1*(2*b(i)+2*(2*cos(theta)^2-1)+2*a(i)*cos(theta));
  C.g(4*i-1) = -0.5*H/K2*(2*c(i)+2*(1+d(i))*cos(theta));
  C.g(4*i)   = -0.5*H/K2*(2*d(i)+2*(2*cos(theta)^2-1)+2*c(i)*cos(theta));
end 
C.g(n) = H/x(n);
C.G = C.g*C.g'/H;
for i = 1:k
  K1 = 1+a(i)^2+b(i)^2+2*b(i)*(2*cos(theta)^2-1)+2*a(i)*(1+b(i))*cos(theta); 
  K2 = 1+c(i)^2+d(i)^2+2*d(i)*(2*cos(theta)^2-1)+2*c(i)*(1+d(i))*cos(theta);
  C.G(4*i-3,4*i-3) = -C.G(4*i-3,4*i-3) + H*(0.5/K1*2);
  C.G(4*i-3,4*i-2) = -C.G(4*i-3,4*i-2) + H*(0.5/K1*2*cos(theta));
  C.G(4*i-2,4*i-3) = C.G(4*i-3,4*i-2);
  C.G(4*i-2,4*i-2) = -C.G(4*i-2,4*i-2) + H*(0.5/K1*2);

  C.G(4*i-1,4*i-1) = 3*C.G(4*i-1,4*i) - H/K2;
  C.G(4*i-1,4*i) = 3*C.G(4*i-1,4*i) - H/K2*cos(theta);
  C.G(4*i,4*i-1) = 3*C.G(4*i,4*i-1) - H/K2*cos(theta);
  C.G(4*i,4*i) = 3*C.G(4*i,4*i-1) - H/K2;
end 
C.G(n,n) = 0;
C.g = flag * C.g;
C.G = flag * C.G;
end 
