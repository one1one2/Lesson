function [x,f,info] = Leastp(fun,Point)
path(path,'linesearch');

x = Point;
n = length(x);

epsilon_f = 1e-8;
epsilon_g = 1e-15;
epsilon = 1e-10;
step = 100;

p = 1;
c = 1.1;
rho = 0.8;
sigma = 0.1;

info.ite=0; 
info.feva=0;
info.err=0;
info.f=zeros(step,1);
info.g=zeros(step,1);

F = fun(x);
m = length(F);
%mu = zeros(m,1);
%mu(1) = 1;
mu = 1/m*ones(m,1);
xi = -50;

[f0 g] = Fun(fun, x, p, xi, mu);
while (1==1)
  if info.ite>=step info.err=1; break; end  
  [f, g, G] = Fun(fun, x, p, xi, mu);
  if norm(g) <= epsilon_g info.err=2; break; end 
  %G = G + 1*eye(size(G,1));
  d = -G\g;
  d=d/norm(d);
  if g'*d > 0   d = -d; end
  %Rule.opt(1) = 0;
  %Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  Rule.opt(1) = 1;
  Rule.opt = bodfltchk(Rule.opt, [1 1e-6 10 rho sigma]);
  %Rule.opt = bodfltchk(Rule.opt, [1 1e-6 10 0.95 0.05]);
  Rule.crtr = @bostwlf;
  Rule.mthd = @bointrplt33;
  [alpha,info1,perf] = bolinesearch(@(x)Fun(fun,x,p,xi,mu),x,d,Rule);
  info.feva=info.feva+info1(3)+1;
  info.ite=info.ite+1;
  f1=perf.F;
   abs(g'*(perf.x-x)) 
  if abs(g'*(perf.x-x)) < epsilon info.err=0; break; end 
  x=perf.x;
  g=perf.g;
  info.f(info.ite) = f1;
  info.g(info.ite) = norm(g);
  %x1 = fminsearch(@(x)Fun(fun,x,pi,xi,mu),x);
  %[BB f1] = fun(x1);
  %f0 - f1
  %x = x1;
  if abs(f0-f1)<=epsilon_f break; end 
  f0=f1
  p = p*c;
  mu = Getu(fun, x, p, xi, mu);
end
f= Fun(fun, x, p, xi, mu);
[AA FF] = fun(x)
end

function [mu] = Getu(fun, x, p, xi, mu)
[F FF]= fun(x);
M = FF - xi;
q = sign(M)*p;
m = length(F);
v = mu;
for i = 1:m
  if (F{i}.f - xi)*M>0
    v(i) = mu(i)*((F{i}.f-xi)/M)^(q-1);
  end 
end 
mu = v/norm(v,1);

end 

function [f, g, G] = Fun(fun, x, p, xi, mu)
[F FF]= fun(x);
M = FF - xi;
q = sign(M)*p;
m = length(F);
n = length(x);
K = 0;
g = zeros(n,1);
G = zeros(n,n);
for i = 1:m
  if (F{i}.f - xi)*M>0
    K = K + mu(i)*(F{i}.f-xi)^q;
  end 
end 
f = K^(1/q);
for j = 1:n
  for i = 1:m
    if (F{i}.f - xi)*M>0
      g(j) = g(j) + K^(1/q-1)*mu(i)*(F{i}.f-xi)^(q-1)*F{i}.g(j);
    end 
  end
end 

if nargout > 2
  for j = 1:n
    for k = 1:n
      for i = 1:m
        if (F{i}.f - xi)*M>0
          G(j,k) = G(j,k) + K^(1/q-1)*mu(i)*(F{i}.f-xi)^(q-1)*F{i}.G(j,k)  ... 
          + K^(1/q-1)*mu(i)*(q-1)*(F{i}.f-xi)^(q-2)*F{i}.g(j)*F{i}.g(k) ...
          + (1/q-1)*K^(1/q-2)*mu(i)^2*(F{i}.f-xi)^(2*q-2)*F{i}.g(j)*F{i}.g(k);
        end 
      end 
    end
  end 
end

end 
