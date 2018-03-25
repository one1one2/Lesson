function [r J] = extend_powell_sglJ(x)
%EXTEND_POWELL_SGLJ fExtended Powell sigular function with output r J
% (3, -1, 0, 1, ..., 3, -1, 0, 1) for test
% min = 0 at (0, ..., 0)

% Version:	2009.05.22
% Create:	2009.05.22
% Coder:	Xin Liang
% Bug Submission:	liangxinslm@163.com

n = length(x);
n = n - mod(n, 4);

for i = 1:4:n-3
    [r(i:i+3,1) J(i:i+3,i:i+3)] = powell_sglJ(x(i:i+3));
end
