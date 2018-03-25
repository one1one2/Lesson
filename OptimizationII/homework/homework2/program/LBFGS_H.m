function [x,f,info] = LBFGS_H(fun,Point)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% 有限内存BFGS方法
% 调用方式：
% [x,f,info]=LBFGS(fun,Point)
% Input: 
%   fun  in P, function_handle,scalar  所求函数
%   Point x in P,n-vector  初始点
% Output:
%   x  x*  所求最小点 in P,n-vector
%   f  f(x*) 最小值  ,scalar
%   info.ite  迭代次数
%   info.feva  函数调用次数
%   info.err  是否满足精度要求  0满足 1不满足
path(path,'linesearch');

x = Point;
n = length(x);
m = 30; 

epsilon_f = 1e-25;
epsilon_g = 1e-15;
step = 10000;

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
  q = g;
  k = min(info.ite, m);
  for i = k:-1:1
    alpha(i) = rho(i)*S(:,i)'*q;
    q = q - alpha(i)*Y(:,i);
  end 
  r = (S(:,k)'*Y(:,k))/(Y(:,k)'*Y(:,k))*q;
  for i = 1:k
    beta = rho(i)*Y(:,i)'*r;
    r = r + S(:,i)*(alpha(i) - beta);
  end 

  d = -r;
else 
  d = -g;
end
  d = d/norm(d);
  %if g'*d > 0   d = -d; end
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
