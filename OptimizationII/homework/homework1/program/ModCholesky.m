function [L D E] = ModCholesky(G)
%% Modified Cholesky Factorization (Gill and Murray)
%% Optimization theory and methods nonlinear programming
%% Page 138 Alg 3.3.2

n = size(G,1);

L = zeros(size(G));
D = zeros(size(G));
E = zeros(size(G));

% step 1
epsilon = 1e-8;
gamma = 0;
xi = 0;
for i = 1:n
  gamma = max(gamma, abs(G(i,i)));
  for j = 1:n
    if (j~=i)
      xi = max(xi, abs(G(i,j)));
    end
  end
end

beta = max(max(gamma, xi/sqrt(n*n-1)), epsilon);
delta = epsilon;
C = G;

% step 2--6
for j = 1:n
  
  % step 2
  %q = j;  
  %maxabs = abs(C(j,j));
  %for i = j:n
    %if (abs(C(i,i)) > maxabs)
      %q = i;
      %maxabs = abs(C(i,i));
    %end
  %end 
  %G([j q],:) = G([q j],:);
  %G(:,[j q]) = G(:,[q j]);

  % step 3
  for s = 1:j-1
    L(j,s) = C(j,s)/D(s,s);
  end 
  L(j,j) = 1;
  theta = 0;
  for i = j+1:n
    C(i,j) = G(i,j);
    for s = 1:j-1
      C(i,j) = C(i,j) - L(j,s)*C(i,s);
    end
    theta = max(theta, abs(C(i,j)));
  end
  % step 4
  D(j,j) = max(max(delta, abs(C(j,j))), theta^2/beta);
  E(j,j) = D(j,j) - C(j,j);
  % step 5 
  C(j,j) = L(j,j)*D(j,j);
  for i = j+1:n
    C(i,i) = C(i,i) - C(i,j)^2/D(j,j);
  end 
end 
end 
