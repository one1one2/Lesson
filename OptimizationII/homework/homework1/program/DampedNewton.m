function [x,f,info] = DampedNewton(fun,Point)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% 带步长牛顿方法
% 调用方式：
% [x,f,info]=DampedNewton(fun,Point)
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
e = 1e-15;
info.ite=0; 
info.feva=0;
info.err=0;
step = 1000;
[f0 g H] = fun(x);
while (1==1) 
  if info.ite>=step info.err=1; break; end
  
  [f, g, H] = fun(x);
  d=-(H\g)';

  Rule.opt(1) = 0;
  Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  %Rule.opt = bodfltchk(Rule.opt, [1 10 10 0.95 0.05]);
  Rule.crtr = @bostwlf;
  Rule.mthd = @bointrplt33;
  [alpha,info1,perf]=bolinesearch(fun,x,d,Rule);
  info.feva=info.feva+info1(3)+1;
  f1=perf.F;
  x=perf.x; 
  g=perf.g;
  if abs(f0-f1)<=e break; end 
  f0=f1;  
  info.ite=info.ite+1; 
end
f=fun(x);
end
