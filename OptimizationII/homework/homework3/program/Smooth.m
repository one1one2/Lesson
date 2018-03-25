function [x,f,info] = Smooth(fun,Point)
path(path,'linesearch');

x = Point;
n = length(x);

epsilon_f = 1e-25;
epsilon_g = 1e-15;
epsilon = 1e-6;
step = 100;

mu = 100;
beta = 0.5;
rho = 0.8;
sigma = 0.1;

info.ite=0; 
info.feva=0;
info.err=0;
info.f=zeros(step,1);
info.g=zeros(step,1);

[f0 g] = Fun(fun, x, mu);
while (1==1)
  if info.ite>=step info.err=1; break; end  
  [f, g, G] = Fun(fun, x, mu);
  if norm(g) <= epsilon_g info.err=2; break; end 
  %G = G + 1e-2*eye(size(G,1));
  d = -G\g;
  %d=d/norm(d);
  if g'*d > 0   d = -d; end
  %Rule.opt(1) = 0;
  %Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  Rule.opt(1) = 1;
  Rule.opt = bodfltchk(Rule.opt, [1 1e-5 100 rho sigma]);
  %Rule.opt = bodfltchk(Rule.opt, [1 1e-6 10 0.95 0.05]);
  Rule.crtr = @bostwlf;
  Rule.mthd = @bointrplt33;
  [alpha,info1,perf] = bolinesearch(@(x)Fun(fun,x,mu),x,d,Rule);
  info.feva=info.feva+info1(3)+1;
  info.ite=info.ite+1;

  f1=perf.F;
   abs(g'*(perf.x-x)) 
  if abs(g'*(perf.x-x)) < epsilon info.err=0; break; end 
  x=perf.x
  g=perf.g;
  info.f(info.ite) = f1;
  info.g(info.ite) = norm(g);
  if abs(f0-f1)<=epsilon_f break; end 
  f0=f1
  mu = mu*beta;
end
f= Fun(fun, x, mu);
end


function [f, g, G] = Fun(fun, x, mu)
[F FF]= fun(x);
m = length(F);
n = length(x);
K = 0;
g = zeros(n,1);
G = zeros(n,n);
for i = 1:m
  K = K + exp((F{i}.f-FF)/mu);
end 
f = FF + mu*log(K);
for j = 1:n
  for i = 1:m
    g(j) = g(j) + exp((F{i}.f-FF)/mu)/K*F{i}.g(j);
  end
end 

if nargout > 2
  for j = 1:n
    for k = 1:n
      for i = 1:m
        G(j,k) = G(j,k) + exp((F{i}.f-FF)/mu)*F{i}.G(j,k)/K  ... 
        + exp((F{i}.f-FF)/mu)*F{i}.g(j)*F{i}.g(k)/mu/K ...
        - exp((F{i}.f-FF)/mu)^2*F{i}.g(j)*F{i}.g(k)/mu/K/K;
      end 
    end
  end 
end

end 
