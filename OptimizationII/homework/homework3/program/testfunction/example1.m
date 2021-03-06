function [A F] = example1(x)
A{1}.f = x(1)^2 + x(2)^4;
A{1}.g = zeros(2,1);
A{1}.g(1) = 2*x(1);
A{1}.g(2) = 4*x(2)^3;
A{1}.G = zeros(2,2);
A{1}.G(1,1) = 2;
A{1}.G(2,2) = 12*x(2)^2;

A{2}.f = (2-x(1))^2  + (2-x(2))^2;
A{2}.g = zeros(2,1);
A{2}.g(1) = 2*(x(1)-2);
A{2}.g(2) = 2*(x(2)-2);
A{2}.G = zeros(2,2);
A{2}.G(1,1) = 2;
A{2}.G(2,2) = 2;

K = 2*exp(x(2)-x(1));
A{3}.f = K;
A{3}.g = zeros(2,1);
A{3}.g(1) = -K;
A{3}.g(2) = K;
A{3}.G = [1 -1; -1 1]*K;

F = max(A{1}.f, A{2}.f);
F = max(F, A{3}.f);
end 


