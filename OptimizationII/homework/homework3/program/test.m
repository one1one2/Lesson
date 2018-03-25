path(path,'testfunction/');

%method = @Smooth;
%method = @Leastp;
method = @SQP;
%flag = 'example3';
flag = 'digital_filter';

if strcmp(flag, 'simpletest')
  x0 = [1,1]';
  [x,f,info] = method(@simpletest,x0)
end 

if strcmp(flag, 'example1')
  x0 = [1,-0.1]';
  [x,f,info] = method(@example1,x0)
end 
if strcmp(flag, 'example2')
  x0 = [1,-0.1]';
  [x,f,info] = method(@example2,x0)
end 
if strcmp(flag, 'example3')
  x0 = [0 0 0 0]';
  [x,f,info] = method(@example3,x0)
end 
if strcmp(flag, 'digital_filter')
  x0 = [0 0.999 0 -0.15 0 -0.68 0 -0.72 0.37]';
  [x,f,info] = method(@digital_filter,x0)
end 
