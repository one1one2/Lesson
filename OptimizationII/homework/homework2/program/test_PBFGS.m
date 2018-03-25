path(path,'testfunction');

method = @PBFGS;
flag = 'extend_wood';


if strcmp(flag, 'LMS')
  n = 30;
  x0 = 9*ones((n-1)^2,1);
  t = cputime;
  [x,f,info] = method(@LMS,n^2,x0)
  cputime - t
end 


if strcmp(flag,'S303')
  n = 100;
  x0 = 0.1*ones(n,1);
  t = cputime;
  [x,f,info] = method(@S303,n+2,x0)
  t = cputime - t
end 

if strcmp(flag,'GenCube')
  n = 10000;
  x0 = 1.5*ones(n,1); 
  t = cputime;
  [x,f,info] = method(@GenCube, n, x0)
  t = cputime - t
end 


if strcmp(flag, 'extend_wood')
  n = 50;
  x0=zeros(n,1);
  for i=1:(n/2)
    x0(2*i-1,1) = -3;
    x0(2*i,1) = -1;
  end
  t = cputime;
  [x,f,info] = method(@extend_wood,n/2-1 ,x0)
  cputime - t
end 

