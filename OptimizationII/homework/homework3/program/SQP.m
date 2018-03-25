function [x,f,info] = SQP(fun,Point)
path(path,'linesearch');

x = Point;
n = length(x);

epsilon_f = 1e-25;
epsilon_g = 1e-15;
epsilon = 1e-8;
step = 100;

delta = 0.5;
beta = 0.8;
sigma = 0.5;

info.ite=0; 
info.feva=0;
info.err=0;

[AA f0] = fun(x);
m = length(AA);
B = eye(n);
while (1==1)
  if info.ite>=step info.err=1; break; end  
  [AA FF] = fun(x);
  H = zeros(n+1);
  H(1:n,1:n) = B;
  H(n+1,n+1) = delta;
  f = zeros(n+1,1);
  f(n+1) = 1;
  A = zeros(m,n+1);
  b = zeros(m,1);
  for i = 1:m 
    A(i,1:n) = AA{i}.g;
    A(i,n+1) = -1;
    b(i) = FF - AA{i}.f;
  end 
  [res,fval,exitflag,output,lambda]=quadprog(H,f,A,b);
  lambda = lambda.ineqlin;
  d = res(1:n);
  t = res(n+1);

  d = d/(1+delta*t);
  lambda = lambda/(1+delta*t);
  if norm(d) < epsilon info.err=0; break; end 
  
  alpha = 1;
  while Fun(fun,x+alpha*d) > FF + sigma*alpha*t
    alpha = alpha*beta;
    info.feva = info.feva + 1;
  end 
  info.ite = info.ite + 1;
  x = x + alpha*d;
  s = alpha*d;
  [BB f1] = fun(x);   
  y = zeros(n,1);
  if abs(f0-f1)<=epsilon_f break; end 
  for i = 1:m
    y = y + lambda(i)*(BB{i}.g - AA{i}.g);
  end
  if y'*s >= 0.2*s'*B*s
    theta = 1;
  else 
    theta = 0.8*s'*B*s/(s'*B*s-y'*s);
  end 
  y = theta*y + (1-theta)*B*s;
  B = B - (B*s)/(s'*B*s)*(s'*B) + y/(y'*s)*y';
  f0=f1
end
f = Fun(fun, x);
end


function [f] = Fun(fun, x)
[A, f]= fun(x);
end 
