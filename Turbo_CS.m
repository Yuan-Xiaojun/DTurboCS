function [nmse,x_hat]=Turbo_CS(y,A,At,params,errx)
% input signal: y
% measurement matrix: A, At
%-------------parameters----------------
M=params.M;
N=params.N;
sigma=params.sigma;
delta=M/N;
iter=params.iter;
lambda = params.lambda;
x_hat=zeros(N,1);
v = norm(y)^2/M-sigma^2;
%xr = params.xr; % real value of x
%-----------iterative algorithm------------------
for ii=1:iter
    r_hat = x_hat+1/delta*At(y-A(x_hat));
    c = (1/delta-1)*v+1/delta*sigma^2;
    % Kernel function to promote sparsity
    [xpost,~,x_hat,v] = mmse_denoiser(r_hat,c,lambda);
    nmse(ii) = errx(xpost); 
end
end