function [x,f,info] = Sorensen(fun,Point)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% More-Sorensen方法 Page160
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

x=Point; 
e_f = 1e-15;
info.ite=0; 
info.feva=0;
info.err=0;
step = 1000;
f0 = fun(x) + 1;
rho = 0.3;
sigma = 0.7;
gamma = 0.7;
while (1==1)
  if info.ite>=step 
    info.err=1; 
    break; 
  end
  [f, g, G] = fun(x);
  info.feva = info.feva+1;

  [s,d] = DescentPair(g,G);
    alpha = 0.99;
    while (1==1)
      x1 = x + alpha^2*s + alpha*d;
      [f1, g1] = fun(x1);
      if (f1 <= f + rho*alpha^2*(s*g + 0.5*d*G*d') && ...
          (d+2*alpha*s)*g1 >= sigma*(d*g + 2*alpha*s*g + alpha*d*G*d'))
        break;
      end
      alpha = alpha*gamma;
      info.feva = info.feva + 1;
    end
  %f1 - f - rho*alpha^2*(s*g + 0.5*d*G*d') 
  %x*g1 - sigma*(d*g + 2*alpha*s*g + alpha*d*G*d')
  x = x1;
  if abs(f0 - f1)<= e_f
    break; 
  end 
  f0 = f1;  
  info.ite = info.ite + 1; 
end
f = fun(x);
end
