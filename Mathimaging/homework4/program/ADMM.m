function u = ADMM(A, f, u_origin)

addpath('./2DTWFT/2D');

frame=3; % type of wavelet frame used: 0 is Haar; 1 is piecewise linear; 3 is piecewise cubic
Level=2; % level of decomposition, typically 1-4.
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel2D(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel2D(x,R,Level); % Frame reconstruction

%设置参数
mu = 0.1;
maxstep = 10;
tol = 1e-7;
delta = 0.1;
tau = 0.1;

%计算一些Fourier变换
A_hat = fft2(A);
f_hat = fft2(f);

Dem = A_hat.*A_hat + mu.*ones(size(A_hat));
Af_hat = A_hat.*f_hat;
f_abs = norm(f, 'fro');

%设置初值

nd = length(D) - 1;
for i = 1:Level
  for j = 1:nd
    for k = 1:nd
      lambda{i}{j,k} = tau*ones(size(f));
    end
  end
end
lambda{1}{1,1} = zeros(size(f));

d = W(f);
b = CoeffOper2D('-', d, d); %b0=0

for iter = 1:maxstep
  % u = (A'A+mu)^{-1}*(A'f+mu*W'(d-b))
  temp = CoeffOper2D('-', d, b);
  u_hat = (Af_hat + mu.*WT(temp))./Dem;
  u = real(ifft2(u_hat));

  % d = shrink(Wu+b)
  temp = CoeffOper2D('+', W(u), b);
  d = CoeffOper2D('vs', temp, lambda);
  %b = b+ delta(Wu-d)
  temp = CoeffOper2D('-', W(u), d); 
  err = cellnorm2D(temp, 2)/f_abs;
  if err < tol
    break;
  end 
  temp = CoeffOper2D('*', temp, delta);
  b = CoeffOper2D('+', b, temp);
end 

err_ADMM = norm(u - u_origin, 'fro')/norm(u_origin ,'fro')
figure
imshow(u,[])




end 
