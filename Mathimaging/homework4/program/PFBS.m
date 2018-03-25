function u = PFBS(A, f, u_origin)


addpath('./2DTWFT/2D');
frame=3; % type of wavelet frame used: 0 is Haar; 1 is piecewise linear; 3 is piecewise cubic
Level=2; % level of decomposition, typically 1-4.
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel2D(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel2D(x,R,Level); % Frame reconstruction

%设置参数
kappa = 1; 
L = 10;%norm(A'*A)+kappa
maxstep = 10;
tol = 1e-8;
tau = 10;

nd = length(D) - 1;
for i = 1:Level
  for j = 1:nd
    for k = 1:nd
      lambda{i}{j,k} = tau*ones(size(f));
    end
  end
end
lambda{1}{1,1} = zeros(size(f));
%设置初值

alpha = W(f);

for iter = 1:maxstep
  temp1 = W(A'*(A*WT(alpha)-f));
  temp2 = CoeffOper2D('-',alpha, W(WT(alpha))); 
  temp3 = CoeffOper2D('*', temp2, kappa);
  temp4 = CoeffOper2D('+', temp1, temp3);
  temp5 = CoeffOper2D('*', temp4, 1/L);
  g = CoeffOper2D('-', alpha, temp5);
  alpha = CoeffOper2D('vs', g, lambda);
end 
u = WT(alpha);
err_ADMM = norm(u - u_origin, 'fro')/norm(u_origin ,'fro')
figure
imshow(u,[])

end 



