function [F g] = LMS(x, flag, F, varargin)
if nargin == 1
   n = sqrt(length(x)) + 1;
   m = n*n;
   F = 0;
   g = zeros(size(x));
   for i = 1:n
     for j = 1:n
       [a b c d] = getABCD(i, j, n, x);
       temp = sqrt(1+m/2*((a-d)^2+(b-c)^2));
       F = F + temp/m;
       if (i>1&&j>1) g((i-2)*(n-1)+j-1) = g((i-2)*(n-1)+j-1) + (a-d)/2/temp; end 
       if (i>1&&j<n) g((i-2)*(n-1)+j) = g((i-2)*(n-1)+j) + (b-c)/2/temp; end 
       if (i<n&&j>1) g((i-1)*(n-1)+j-1) = g((i-1)*(n-1)+j-1) + (c-b)/2/temp; end 
       if (i<n&&j<n) g((i-1)*(n-1)+j) = g((i-1)*(n-1)+j) + (d-a)/2/temp; end 
     end 
   end 
end 
if nargin==2
   n = sqrt(length(x)) + 1;
   m = n*n;
   for i = 1:n
     for j = 1:n
       F{(i-1)*n+j} = zeros(2,1);
       [a b c d] = getABCD(i, j, n, x);
       temp = sqrt(1+m/2*((a-d)^2+(b-c)^2));
       if (i>1&&j>1) 
         F{(i-1)*n+j}(1) = a; 
       end 
       if (i<n&&j<n) 
         F{(i-1)*n+j}(1) = F{(i-1)*n+j}(1) - d;
       end 
       if (i>1&&j<n) 
         F{(i-1)*n+j}(2) = b;
       end 
       if (i<n&&j>1) 
         F{(i-1)*n+j}(2) = F{(i-1)*n+j}(2) - c;
       end 
       g{(i-1)*n+j} = [a-d;b-c]/2/temp;
     end 
   end 
end 
if nargin==3
  n = sqrt(length(F)) + 1;
  for i = 1:n
    for j = 1:n
       if (i>1&&j>1) F((i-2)*(n-1)+j-1) = F((i-2)*(n-1)+j-1) ...
         + x{(i-1)*n+j}(1); end 
       if (i>1&&j<n) F((i-2)*(n-1)+j) = F((i-2)*(n-1)+j) ...
         + x{(i-1)*n+j}(2); end 
       if (i<n&&j>1) F((i-1)*(n-1)+j-1) = F((i-1)*(n-1)+j-1) ...
         - x{(i-1)*n+j}(2); end 
       if (i<n&&j<n) F((i-1)*(n-1)+j) = F((i-1)*(n-1)+j) ...
         - x{(i-1)*n+j}(1); end  
    end 
  end 
end 


end 

function [a b c d] = getABCD(i, j, n, x)
  if (i==1||j==1) 
    a = boundary((j-1)/n,(i-1)/n);
  else 
    a = x((i-2)*(n-1)+j-1);
  end
  if (i==1||j==n)
    b = boundary(j/n,(i-1)/n);
  else
    b = x((i-2)*(n-1)+j);
  end
  if (i==n||j==1)
    c = boundary((j-1)/n,i/n);
  else   
    c = x((i-1)*(n-1)+j-1);
  end 
  if (i==n||j==n)
    d = boundary(j/n,i/n);
  else 
    d = x((i-1)*(n-1)+j);
  end 
end 

function f = boundary(x, y)
f = 4*x - 8*y + 9;
end 
