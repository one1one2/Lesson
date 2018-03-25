% function to get Hessian matrix at point x
% x is a row vector with size(x) = [1,n];

function G = Hesse(fun, x)
n = size(x,2);
G = zeros(n,n);
h = 1e-5;
for i = 1:n
  x1 = x;
  x2 = x;
  x1(i) = x1(i) + h;
  x2(i) = x2(i) - h;
  G(i,i) = ((fun(x1) - fun(x))/h + (fun(x2) - fun(x))/h)/h;
  for j = 1:n
    if (j~=i)
      xne = x1; 
      xnw = x1;
      xne(j) = x1(j) + h;
      xnw(j) = x1(j) - h;
      xse = 2.*x - xnw;
      xsw = 2.*x - xne;
      G(i,j) = ((fun(xne) - fun(xnw))/h - (fun(xse) - fun(xsw))/h)/h;
    end
  end
end
      
