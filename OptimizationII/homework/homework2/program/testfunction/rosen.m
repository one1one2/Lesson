function [F g G] = rosen(x)

[r J] = rosenJ(x);
F = r'*r;
g = 2*J'*r;
if nargout > 2
    G = J'*J;
    G(1,1) = G(1,1) - 20 * r(1);
    G = G*2;
end
