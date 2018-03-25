function [F g G] = extend_powell_sgl(x)
%EXTEND_POWELL_SGL Extended Powell sigular function with output F g G
% 
% (3, -1, 0, 1, ..., 3, -1, 0, 1) for test
% min = 0 at (0, ..., 0)

% Version:	2009.05.22
% Create:	2009.05.22
% Coder:	Xin Liang
% Bug Submission:	liangxinslm@163.com

[r J] = extend_powell_sglJ(x);
F = r'*r;
g = 2*J'*r;
if nargout > 2
    n = length(x);
    n = n - mod(n, 4);
    G = J'*J;
    for i = 1:4:n-3
        G(i+1:i+2,i+1:i+2) = G(i+1:i+2,i+1:i+2) + [2 -4; -4 8]*r(i+2);
        G(i:3:i+3,i:3:i+3) = G(i:3:i+3,i:3:i+3) + ...
            [1 -1; -1 1]*2*sqrt(10)*r(i+3);
    end
    G = G*2;
end
