function [F g] = S303(x, flag, F, varargin)
  n = length(x);
  if nargin==1
    F = 0;
    s = 0;
    g = 2*x;

    for i = 1:n
      s = s+ i*x(i)/2;
      F = F + x(i)^2;
    end 
    F = F + s^2 + s^4;
    for i = 1:n
      g(i) = g(i) + s*i + 2*s^3*i;
    end
  end
  if nargin==2
    for i = 1:n
      F{i} = [x(i)];
      g{i} = 2*[x(i)];
    end 
    summ = 0;
    for j = 1:n  
      summ = summ + j*x(j)/2;
    end 
    F{n+1} = [summ];
    g{n+1} = 2*[summ];
    F{n+2} = [summ];
    g{n+2} = 4*[summ^3];
  end 
  if nargin==3
    n = length(F);
    for i = 1:n
      F(i) = F(i) + x{i};
    end 
    for j = 1:n
      F(j) = F(j) + j/2*x{n+1};
      F(j) = F(j) + j/2*x{n+2};
    end 
  end 
end 

