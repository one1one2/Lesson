path(path,'testfunction/');

method = @LBFGS_H;
%method = @LBFGS_B;
%method = @LSR1;
flag = 'extend_wood';

if strcmp(flag, 'LMS')
  n = 30;
  x0 = 9*ones((n-1)^2,1);
  t = cputime;
  [x,f,info] = method(@LMS,x0)
  cputime - t
end 


if strcmp(flag, 'simpletest')
  x0 = [1,1]';
  [x,f,info] = method(@simpletest,x0)
end 

if strcmp(flag,'GenSinev')
  n = 10000;
  x0 = -1*ones(n,1);
  x0(1) = 4.712389;
  [x,f,info] = method(@GenSinev,x0)
end 

if strcmp(flag,'extend_FreRoth')
  n = 10000;
  x0 = -2*ones(n,1);
  [x,f,info] = method(@extend_FreRoth,x0)
end 

if strcmp(flag,'S303')
  n = 1000;
  x0 = 0.1*ones(n,1);
  t = cputime;
  [x,f,info] = method(@S303,x0)
  t = cputime - t
end 

if strcmp(flag,'GenCube')
  n = 10000;
  x0 = 1.5*ones(n,1); 
  t = cputime;
  [x,f,info] = method(@GenCube,x0)
  t = cputime - t
end 

if strcmp(flag, 'GenPowsg')
  n = 10000;
  x0 = zeros(n,1);
  for i = 1:n/2-1
    x0(2*i-1) = 3;
    x0(2*i) = -1;
  end 
  [x,f,info] = method(@GenPowsg,x0)
end 


if strcmp(flag, 'extend_wood')
  n = 50;
  x0=zeros(n,1);
  for i=1:(n/2)
    x0(2*i-1,1) = -3;
    x0(2*i,1) = -1;
  end
  t=cputime; 
  [x,f,info] = method(@extend_wood,x0)
  cputime - t
end 


if strcmp(flag, 'Biggs')
  x0 = [1,2,1,1,1,1]';
  [x,f,info] = method(@Biggs,x0)
end 

if strcmp(flag, 'rosen')
  x0 = [-1.2,1]';
  [x,f,info] = method(@rosen,x0)
end 

if strcmp(flag, 'watson')
  n = 18;
  x0=zeros(1,n);
  [x,f,info] = method(@watson,x0)
end 


if (strcmp(flag,'dis_boundval'))
  n = 30;
  x0 = zeros(n,1);
  for i = 1:n
    x0(i) = i/(n+1)*(i-n-1)/(n+1);
  end 
  [x,f,info]=method(@dis_boundval,x0)
end 

if (strcmp(flag,'extend_powell_sgl'))
  n = 20;
  x0=zeros(n,1);  
  for i=1:4:n-3
    x0(i)=3;
    x0(i+1)=-1;
    x0(i+2)=0;
    x0(i+3)=1;
  end 
  [x,f,info]=method(@extend_powell_sgl,x0)
end 
