function [r J] = watsonJ(x)
%WATSONJ Watson function with output r J
% 
% 2 <= n <= 31
% (0,..., 0) for test
% min = 2.28767e-3 if n = 6
% min = 1.39976e-6 if n = 9
% min = 4.72238e-10 if n = 12

% TODO %    not finished
% Version:	2009.05.22
% Create:	2009.05.22
% Coder:	Xin Liang
% Bug Submission:	liangxinslm@163.com

%error('not designed now');
n = length(x);
if n < 2
    error('n should be in [2, 31]');
elseif n > 31
    x = x(1:31); 
end

r = zeros(31,1);
J = zeros(31,n);
for i = 1:29 
  sum1 = 0;
  sum2 = 0;
  for j = 1:n
    sum1 = sum1 + (j-1)*x(j)*(i/29)^(j-2);
    sum2 = sum2 + x(j)*(i/29)^(j-1);
  end 
  r(i) = sum1 - sum2*sum2 - 1;
  for j = 1:n
    J(i,j) = (j-1)*(i/29)^(j-2) - 2*sum2*(i/29)^(j-1);
  end 
end 
r(30) = x(1);
r(31) = x(2) - x(1)^2 - 1;
J(30,1) = 1;
J(31,1) = -2*x(1);
J(31,2) = 1;


