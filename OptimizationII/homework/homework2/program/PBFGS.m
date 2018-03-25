function [x,f,info] = PBFGS(fun,m,Point)
path(path,'linesearch');

x = Point;
n = length(x);
% num 为每个分块对应的秩，共m个分块

epsilon_f = 1e-25;
epsilon_g = 1e-15;
step = 1000;

info.ite=0; 
info.feva=0;
info.err=0;

[s, y] = fun(x,0);
for i = 1:m
  B{i} = eye(length(s{i}));
end 

[f0 g] = fun(x);
while (1==1)
  if info.ite>=step info.err=1; break; end  
  [f, g] = fun(x);
  if norm(g) <= epsilon_g info.err=2; break; end
  if (info.ite > 0)
    %norm(Bfun(zeros(size(x)),m,fun,B))
    d = -cgs(@(vec)Bfun(vec,m,fun,B),g,1e-6,10000);
  else 
    d = -g;
  end
  d = d/norm(d);
  if g'*d > 0   d = -d; end

  [alpha, info1, perf] = getStepsize(fun, x, d);

  info.feva=info.feva+info1(3)+1;
  info.ite=info.ite+1; 
  x0 = x;
  f1=perf.F; x=perf.x;  g=perf.g;
  info.f(info.ite) = f1;
  info.g(info.ite) = norm(g);
  if abs(f0-f1)<=epsilon_f break; end 
  f0=f1
  B = updateB(B, fun, m, x0, x);
end
f=fun(x);
end

function [alpha, info, perf] = getStepsize(fun, x, d)
  %Rule.opt(1) = 0;
  %Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  Rule.opt(1) = 1;
  Rule.opt = bodfltchk(Rule.opt, [1 1e-4 10 0.95 0.05]);
  Rule.crtr = @boarmgld;
  Rule.mthd = @bointrplt33;
  [alpha,info,perf] = bolinesearch(fun,x,d,Rule);
end 

function B = updateB(B, fun, m, x0, x)
  [s1 y1] = fun(x0,0);
  [s2 y2] = fun(x,0);
  for i = 1:m
    s = s2{i} - s1{i};
    y = y2{i} - y1{i};
    temp = B{i}*s;
    if norm(temp)>1e-10
      B{i} = B{i} - (temp*temp')/(s'*temp);
    end 
    if norm(y)>1e-10 
      B{i} = B{i} + y*y'/(y'*s);
    end
  end
end 

function res = Bfun(s,m,fun,B)
  res = zeros(size(s)); 
  temp = fun(s,0);
  for i = 1:m
    temp{i} = B{i}*temp{i}; 
  end 
  res = fun(temp, 0, res);
  res = res + s*1e-8;
end 

