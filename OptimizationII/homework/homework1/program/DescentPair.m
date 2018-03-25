function [s d] = DescentPair(g, G)
% 生成下降对，Page162-163

e_g = 1e-8; 
n = size(g,1);
e = 1e-15;
[L, D] = ldl(G);

[V, Lam] = eig(D);

minlambda = min(diag(Lam));
[v, l] = eigs(D,1,'sa');
maximum = max(abs(diag(Lam)));
for i = 1:n
  Lam(i,i) = max(max(abs(Lam(i,i)), e*n*maximum), e);
end
barD = V*Lam*V';
s = -(inv(L')*inv(barD)*inv(L)*g)';
if (minlambda<0 && norm(s) < e_g)
  d = (inv(L')*v)';
else
  d = zeros(1,n);
end
%if (d*g > 0)
  %d = -d;
%end 
%if (s*g > 0)
  %s = -s;
%end 
end 
