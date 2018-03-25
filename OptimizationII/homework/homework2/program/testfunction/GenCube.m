function [F g] = GenCube(x, flag, F, varargin)
n = length(x);
if nargin==1
g = zeros(n,1);
F = (x(1) - 1)^2;
g(1) = 2*(x(1)-1);
for i = 2:n
  F = F + 100*(x(i)-x(i-1)^3)^2; 
  g(i) = g(i) + 200*(x(i)-x(i-1)^3);
  g(i-1) = g(i-1) - 200*(x(i)-x(i-1)^3)*(3*x(i-1)^2); 
end 
end 
%% U 
if nargin == 2
    F{1} = [x(1)];
    g{1} = [2*x(1)-1];
  for i = 2:n
    F{i} = [x(i-1);x(i)];
    g{i} = 200*(x(i) - x(i-1)^3)*[-3*x(i-1)^2;1];
  end 
end 
%% UT
if nargin == 3 
  %F = r;
  F(1) = F(1) + x{1};
  for i = 2:n
    F(i-1:i) = F(i-1:i) + x{i};
  end 

end 

end 


