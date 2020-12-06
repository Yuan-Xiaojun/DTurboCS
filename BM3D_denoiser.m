%denoiser using BM3D
%you should install bm3d first
function [F,F_div] = BM3D_denoiser(y,v)
[m,n]=size(y);

n=sqrt(m*n);
denoiser='fast-BM3D'; %denoiser='BM3D';
denoi=@(noisy,sigma_hat) denoise(noisy,sigma_hat,n,n,denoiser);

denoised=denoi(y,v);
div=0;
K=1;
for ii=1:K
    epsilon=max(y)/100+eps; %epsilon = 0.05*min(v, mean(abs(y),1)) + eps;
    eta=randn(1,n^2); %eta = sign(randn(1,n^2)); % random direction
    rhat_perturbed = y + epsilon*eta';
    xhat_perturbed=denoi(rhat_perturbed,v);
    divp=eta*(xhat_perturbed-denoised)/epsilon;
    div=div+divp;
end
div=div/K/n^2;

F=reshape(denoised,1,n^2);
F_div=div;
end