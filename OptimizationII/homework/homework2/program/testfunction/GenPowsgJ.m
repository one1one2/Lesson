function [r J] = GenPowsgJ(x)

n = length(x);
r = zeros(2*n-4,1);
J = zeros(2*n-4,n);

for i = 1:n/2-1
  [r(4*i-3:4*i) J(4*i-3:4*i,2*i-1:2*i+2)] =...
  powell_sglJ(x(2*i-1:2*i+2)); 
end 
