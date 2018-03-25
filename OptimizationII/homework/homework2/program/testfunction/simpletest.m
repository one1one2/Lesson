function [f, g, H] = simpletest(x)

  f = x(1)^2 + x(2)^2;
  g = 2*x;
  H = 2*eye(2);

end 
