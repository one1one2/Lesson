function [x,f,info] = stableNewton(fun,Point)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% 稳定牛顿方法 Page 150  Algorithm 3.5.4
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
e_f = 1e-15;
e_g = 1e-8;
info.ite = 0; 
info.feva = 0;
info.err = 0;
step = 1000;
f0 = fun(x) + 1;
while (1==1)
  if info.ite>=step 
    info.err=1; 
    break; 
  end
  % step 2
  [f, g, G] = fun(x);
  % step 3
  [L, D, E] = ModCholesky(G);

  if norm(g) > e_g
    % step 4
    d = -(inv(L')*inv(D)*inv(L)*g)';
  else 
    % step 5 Alg3.5.2 
    minimum = D(1,1) - E(1,1);
    flag = 1;
    n = size(D,1);
    for j = 2:n
      if minimum > D(j,j) - E(j,j)
        flag = j;
        minimum = D(j,j) - E(j,j);
      end 
    end
    if (minimum >= 0) 
      break; 
    end 
    e_t = zeros(size(L,1),1);
    e_t(flag) = 1;
    d = (L'\e_t)';
    if (d*g > 0)
      d = -d;
    end 
  end
  Rule.opt(1) = 0;
  Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
  %Rule.opt = bodfltchk(Rule.opt, [1 10 10 0.95 0.05]);
  Rule.crtr = @bostwlf;
  Rule.mthd = @bointrplt33;
  % step 6
  [alpha,info1,perf] = bolinesearch(fun, x, d,Rule);
  info.feva = info.feva + info1(3) + 1;
  f1 = perf.F;   
  x = perf.x; 
  g = perf.g;
  if abs(f0 - f1)<= e_f
    break; 
  end 
  f0 = f1;  
  info.ite = info.ite + 1; 
end
f = fun(x);
end
