function A = simpletest(x)
A{1}.f = x'*x;
A{1}.g = 2*x;
A{1}.G = 2*eye(length(x));
A{2}.f = x'*x/2 + x(1);
A{2}.g = x;
A{2}.g(1) = A{2}.g(1) + 1;
A{2}.G = eye(length(x));
end 
