% function to get gradient at x 
% x is a row vector with size(x) = [1,n];
function g = Gradient(fun, x)
  h = 1e-10;
  n = size(x,2);
  for i = 1:n 
    x1 = x;
    x2 = x;
    x1(i) = x(i) + h;
    x2(i) = x(i) - h;
    g(i) = (fun(x1) - fun(x2))/2/h;
  end 
end 
