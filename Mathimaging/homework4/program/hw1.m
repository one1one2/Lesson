clear
sigma1 = 1.5;

%原图
u_origin = rgb2gray(im2double(imread('Lenna.jpg')));
sigma2 = max(max(u_origin))/10;

%生成待处理图像
ker = fspecial('gaussian',[15 15],sigma1);
tmp = zeros(size(u_origin));
tmp(1:size(ker,1),1:size(ker,2)) = ker;
A = circshift(tmp,-floor(size(ker)/2));
%加模糊，噪声
f = real(ifft2(fft2(u_origin).*fft2(A)));
f = f + sigma2 * randn(size(f));

figure 
imshow(f,[])
%参数

mu = 1;
lambda = 0.001;
maxstep = 100;
tol = 1e-7;
delta = 0.1;

%求解

Dx_ker = [0,0,0;-1./2,0,1./2;0,0,0];
Dy_ker = [0,1./2,0;0,0,0;0,-1./2,0];
L_ker = [0,1,0;1,-4,1;0,1,0];
L = zeros(size(f));
L(1:3,1:3) = L_ker;
L = circshift(L,[-1,-1]);
%定义辅助函数，求导
Dx = @(f) imfilter(f, Dx_ker, 'circ');
Dy = @(f) imfilter(f, Dy_ker, 'circ');

%计算一些Fourier变换
A_hat = fft2(A);
f_hat = fft2(f);
L_hat = fft2(L);

D = A_hat.*A_hat - mu*L_hat;
Af_hat = A_hat.*f_hat;
f_abs = norm(f);

step = 0;
err = 1;
dx = Dx(f);
dy = Dy(f);
bx = dx; by = dy;
% 初始误差
err_true = norm(f - u_origin)/norm(u_origin)
%开始迭代
while err > tol && step < maxstep
	step = step + 1;
	u_hat = (Af_hat - mu*fft2(Dx(dx-bx)+Dy(dy-by)))./D;
	u = real(ifft2(u_hat));

	Dxu = Dx(u); Dyu = Dy(u);
	tempx = Dxu + bx; tempy = Dyu + by;
    para = abs(tempx) + abs(tempy);
    dx = tempx./para.*(max(para - lambda/mu, 0));
	dy = tempy./para.*(max(para - lambda/mu, 0));
	bx = bx + delta*(Dxu - dx);
	by = by + delta*(Dyu - dy);

	errx = norm(Dxu - dx)/f_abs;
	erry = norm(Dyu - dy)/f_abs;
    err = sqrt(errx^2 + erry^2);
end
figure
imshow(u, []);
err = norm(u - u_origin)/norm(u_origin)



