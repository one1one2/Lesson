function [r J] = BiggsJ(x)

m = 13;

r = zeros(m,1);
J = zeros(m,6);
for i = 1:m
  t = 0.1*i;
  r(i) = x(3)*exp(-t*x(1))-x(4)*exp(-t*x(2))...
        +x(6)*exp(-t*x(5))-exp(-t)+5*exp(-10*t)-3*exp(-4*t);
  J(i,1) = x(3)*exp(-t*x(1))*(-t);
  J(i,2) = -x(4)*exp(-t*x(2))*(-t);
  J(i,3) = exp(-t*x(1));
  J(i,4) = -exp(-t*x(2));
  J(i,5) = x(6)*exp(-t*x(5))*(-t);
  J(i,6) = exp(-t*x(5));
end 
end 
