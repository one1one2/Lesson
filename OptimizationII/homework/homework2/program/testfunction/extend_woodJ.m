function [r J] = extend_woodJ(x)

n = length(x);

for i = 1:(n-2)/2
    [r(6*i-5:6*i,1) J(6*i-5:6*i,2*i-1:2*i+2)] = woodJ(x(2*i-1:2*i+2));
end
