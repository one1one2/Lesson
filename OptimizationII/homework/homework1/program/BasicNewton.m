function [x,f,info] = BasicNewton( fun, Point)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% 基本牛顿方法
% 调用方式：
% [x,f,info]=BasicNewton(fun,Point)
% Input: 
%   fun  in P, function_handle,scalar  所求函数
%   Point x in P,n-vector  初始点
% Output:
%   x  x*  所求最小点 in P,n-vector
%   f  f(x*) 最小值  ,scalar
%   info.ite  迭代次数
%   info.feva  函数调用次数
%   info.err  是否满足精度要求  0满足 1不满足
x=Point; 
e=1e-8;
f0=fun(x)+1;
info.ite=0; 
info.feva=0;
info.err=0;
step=1000;
g=g_func(fun,x);
while (norm(g)>e) 
	if info.ite>=step info.err=1; break; end
    H=H_func(fun,x);
    d=-(inv(H)*g)';
    info.feva=info.feva+1;
	x=x+d;
	f1=fun(x);
	g=g_func(fun,x);
	if abs(f0-f1)<=e break; end 
    f0=f1;  info.ite=info.ite+1; 
end
f=fun(x);
end
