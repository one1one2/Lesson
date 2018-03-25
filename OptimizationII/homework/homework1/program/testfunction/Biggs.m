function [F g G] = Biggs(x)

[r J] = BiggsJ(x);
m = size(r,1);
F = r'*r;
g = 2*J'*r;
if nargout > 2
    G = J'*J;
    for i = 1:m 
      t = 0.1*i;
      G(1,1) = G(1,1) + t*t*x(3)*exp(-t*x(1))*r(i);
      G(2,2) = G(2,2) - t*t*x(4)*exp(-t*x(2))*r(i);
      G(5,5) = G(5,5) + t*t*x(6)*exp(-t*x(5))*r(i);
      G(1,3) = G(1,3) - t*exp(-t*x(1))*r(i);
      G(2,4) = G(2,4) + t*exp(-t*x(2))*r(i);
      G(5,6) = G(5,6) - t*exp(-t*x(5))*r(i);
      G(3,1) = G(1,3);
      G(4,2) = G(2,4);
      G(6,5) = G(5,6);
    end 
    G = G*2;
end
