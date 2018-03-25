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

imshow(u_origin,[])
figure 
imshow(f,[])
%参数

err_init = norm(f-u_origin, 'fro')/norm(u_origin, 'fro')
%TV(A, f, u_origin);
%ADMM(A, f, u_origin);
PFBS(A, f, u_origin);
