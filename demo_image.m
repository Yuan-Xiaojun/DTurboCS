%----- D-Turbo-CS Compressive Image Sensing Demo ------
%- paper: Denoising-Based Turbo Compressed Sensing -
%---- Zhipeng Xue, Junjie Ma, Xiaojun Yuan ---------
%---------------------------------------------------
clear all;
close all;
warning('off');
%------------------load images----------------------
img = imread('lena.bmp');
% img = imread('Catch.jpg');
% img = imread('boat.png');
% img = imread('barbara.png');
% img = imread('fingerprint.png');
% img = imread('house.png');
% img = imread('camera.bmp');
% img = imread('barbara.png');
% img = imread('brain.png');
% img = imread('phantom.png');
[D,~] = size(img); % square image
img = double(img); %img = double(rgb2gray(img));
%------------------parameters------------------------
rate = 0.05; % measurement rate
M = floor(D^2*rate); % number of measurment
type = "randsecdctsign"; % measurement type
[A,At] = LinerOperator(D^2,M,type);
%-----------image measurement------------------------
sigma = 0;
x_0 = reshape(img,D^2,1);
x_0 = dct(x_0); % if K is choosen as 1, this line should be deleted
y = A(x_0)+sigma*randn(M,1);
errx = @(z)norm(z-x_0)^2/norm(x_0)^2;
%--------------reconstruction------------------------
params.M = M;
params.N = D^2;
params.sigma = sigma;
% denoiser type: set to 1 only if BM3D denoiser is installed, for sure-let
% denoider set to 3
params.K = 3;
[errx_out,x_hat] = TurboCS_img(y,A,At,params,errx);
x_out= idct(x_hat);a
img_recon = reshape(x_out,D,D);
%----------------------figures------------------------
figure(1);
imshow(uint8(img_recon));
Img_Max=255;
MSE = 1/D^2*sum(sum((img - img_recon).^2));
PSNR = 10*log10(Img_Max^2/MSE);
disp(['Recovery PSNR : ',num2str(PSNR)]);