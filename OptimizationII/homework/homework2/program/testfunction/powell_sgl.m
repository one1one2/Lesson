function [F g G] = powell_sgl(x)

[r J] = powell_sglJ(x);
F = r'*r;
g = 2*J'*r;
if nargout > 2
    G = J'*J;
    G(1,1) = G(1,1) + 2*sqrt(10)*r(4);
    G(2,2) = G(2,2) + 2*r(3);
    G(3,3) = G(3,3) + 8*r(3);
    G(4,4) = G(4,4) + 2*sqrt(10)*r(4);
    G = G*2;
end
