clear

%原图
u_origin = imresize(rgb2gray(im2double(imread('baby.jpg'))),1);
%u_origin = rgb2gray(im2double(imread('Lenna.jpg')));

%生成待处理图像
sigma1 = 5;
ker = fspecial('gaussian',[15 15],sigma1);
tmp = zeros(size(u_origin));
tmp(1:size(ker,1),1:size(ker,2)) = ker;
A = circshift(tmp,-floor(size(ker)/2));
%加模糊，噪声
u_origin = real(ifft2(fft2(u_origin).*fft2(A)));

%加噪声
%sigma2 = max(max(u_origin))/10;
%f = u_origin + sigma2 * randn(size(u_origin));

imshow(u_origin,[]);
figure

%参数
dt = 0.5;
step = 300; 
alpha = 1;
tol = 1e-7;
Reinit_step = 10;

u = levelset(u_origin, dt, step, alpha, tol, Reinit_step);

