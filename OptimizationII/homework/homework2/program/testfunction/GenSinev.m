function [F g] = GenSinev(x)
n = length(x);
F = 0;
c1 = 1e-4;
c2 = 4;
g = zeros(n,1);
for i = 1:n-1
  F = F + (x(i+1)-sin(x(i)))^2/c1 + x(i)^2/c2;
  g(i) = g(i) + 2/c1*(x(i+1)-sin(x(i)))*(-cos(x(i))) + 2*x(i)/c2;
  g(i+1) = g(i+1) + 2/c1*(x(i+1)-sin(x(i)));
end 

end 
