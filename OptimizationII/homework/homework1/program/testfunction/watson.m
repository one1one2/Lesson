function [F g G] = watson(x)
%WATSON Watson function with output F g G
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
[r J] = watsonJ(x);
n = length(x);
F = r'*r;
g = 2*J'*r;
if nargout > 2
    G = J'*J;
    for i = 1:29
      for j = 1:n
        for k = 1:n
          G(j,k) = G(j,k) - 2*(i/29)^(j+k-2)*r(i);
        end 
      end
    end 
    G(1,1) = G(1,1) - 2*r(31);
	G = G*2;
end
