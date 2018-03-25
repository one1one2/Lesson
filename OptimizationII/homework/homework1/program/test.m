path(path,'testfunction/');

%method = @DampedNewton;
method = @StableNewton;
%method = @Sorensen;
flag = 'watson';

if strcmp(flag, 'Biggs')
  x0 = [1,2,1,1,1,1];
  [x,f,info] = method(@Biggs,x0)
end 

if strcmp(flag, 'rosen')
  x0 = [-1.2,1];
  [x,f,info] = method(@rosen,x0)
end 

if strcmp(flag, 'watson')
  n = 18;
  x0=zeros(1,n);
  [x,f,info] = method(@watson,x0)
end 


if (strcmp(flag,'dis_boundval'))
  n = 30;
  x0 = zeros(1,n);
  for i = 1:n
    x0(i) = i/(n+1)*(i-n-1)/(n+1);
  end 
  [x,f,info]=method(@dis_boundval,x0)
end 

if (strcmp(flag,'extend_powell_sgl'))
  x0=zeros(1,n);  
  for i=1:4:n-3
    x0(i)=3;
    x0(i+1)=-1;
    x0(i+2)=0;
    x0(i+3)=1;
  end 
  [x,f,info]=method(@extend_powell_sgl,x0)
end 
