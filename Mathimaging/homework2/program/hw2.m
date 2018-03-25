clear

%原图
u_origin = rgb2gray(im2double(imread('Lenna.jpg')));
sigma2 = max(max(u_origin))/10;


%加噪声
f = u_origin + sigma2 * randn(size(u_origin));


figure 
imshow(u_origin,[]);
figure
imshow(f,[]);

%参数
t_end = 1;
dt = 0.01;

step = t_end/dt;

%热方程 
%u_heat = heat(f, dt, step, u_origin, 1);
%err = norm(f - u_origin)/norm(u_origin)
%figure 
%imshow(u_heat, []);
%err_heat = norm(u_heat - u_origin)/norm(u_origin)

%PM方程
%u_pm = pm(f, dt, step, u_origin);
%figure 
%imshow(u_pm, []);
%err_pm = norm(u_pm - u_origin)/norm(u_origin)

%SF方程
%用热方程来生成模糊图片，再用sf方法去模糊
u_heat = heat(f, dt, step, u_origin, 1);
figure 
imshow(u_heat, []);
%
u_shock = shock_filter(u_heat, 10*dt, 3*step, u_origin);
figure 
imshow(u_shock, []);
err_shock = norm(u_shock - u_origin)/norm(u_origin)

