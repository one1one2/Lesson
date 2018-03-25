function [x,f,info] = LSR1(fun,Point)
path(path,'linesearch');

x = Point;
n = length(x);
m = 3; 

epsilon_f = 1e-25;
epsilon_g = 1e-15;
step = 100000;

info.ite=0; 
info.feva=0;
info.err=0;
info.f=zeros(step,1);
info.g=zeros(step,1);

S = zeros(n, m);
Y = zeros(n, m);
rho = zeros(1, m);
alpha = zeros(1, m);

[f0 g] = fun(x);
while (1==1)
  if info.ite>=step info.err=1; break; end  
  [f, g] = fun(x);
  if norm(g) <= epsilon_g info.err=2; break; end 
  
if (info.ite > 0)
  k = min(info.ite, m);
  gamma = (S(:,k)'*Y(:,k))/(Y(:,k)'*Y(:,k));
  %gamma = 100;
  r = gamma*g;
  temp = S(:,1:k)'*g-Y(:,1:k)'*g*gamma;
  temp = (1e-8*eye(k)+(S(:,1:k)-Y(:,1:k)*gamma)'*Y(:,1:k))\temp;
  temp = S(:,1:k)*temp-Y(:,1:k)*temp*gamma;
  d = - r - temp;
else 
  d = -g;
end
  d = d/norm(d);
  if g'*d > 0   d = -d; end
  %Rule.opt(1) = 0;
  %Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  Rule.opt(1) = 1;
  Rule.opt = bodfltchk(Rule.opt, [1 1e-4 10 0.95 0.05]);
  Rule.crtr = @boarmgld;
  Rule.mthd = @bointrplt33;
  [alpha,info1,perf] = bolinesearch(fun,x,d,Rule);
  info.feva=info.feva+info1(3)+1;

  info.ite=info.ite+1; 
  s = perf.x - x;
  y = perf.g - g;
  if info.ite <= m 
    S(:,info.ite) = s;
    Y(:,info.ite) = y;
    rho(info.ite) = 1/(s'*y);
  else 
    S(:,1:m-1) = S(:,2:m);
    S(:,m) = s;
    Y(:,1:m-1) = Y(:,2:m);
    Y(:,m) = y;
    rho(:,1:m-1) = rho(:,2:m);
    rho(m) = 1/(s'*y);
  end 
  

  f1=perf.F;
  x=perf.x;
  g=perf.g;
  info.f(info.ite) = f1;
  info.g(info.ite) = norm(g);
  if abs(f0-f1)<=epsilon_f break; end 
  f0=f1
end
f=fun(x);
end
