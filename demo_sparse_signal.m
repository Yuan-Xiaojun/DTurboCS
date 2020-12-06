% denoising based turbo-cs

clear all;
clc;

% signal and measurement parameters
rou=0.27;
N=20000;
delta=0.5;
iter=70;
M=delta*N;
M=fix(M);

% generation of sparse signal x
pattern=random('Binomial',1,rou,N,1);
x = 1/sqrt(rou)*randn(N,1);
x = x.* pattern;

% linear measurement operator
[A, At] = LinerOperator(N,M,'highrandsecdct');
% [A, At] = LinerOperator(N,M,'gaurand');

% measurement
sigma=sqrt(M/N)*sqrt((10^(-50/10)));
y = A(x)+sigma*randn(M,1);

% normlized mse
errx = @(z)norm(z-x)^2/N;

% algorithm parameters
params.M = M;
params.N = N;
params.sigma = sigma;
params.iter = iter;
params.K = 3;
params.xr = x;
params.lambda = rou;
params.type = "orthogonal";

% dturbocs algorithm
[errx_out1,x_hat1]=DTurboCS(y,A,At,params,errx);
[errx_out2,x_hat2]=SURE_AMP(y,A,At,params,errx);
[errx_out3,x_hat3]=Turbo_CS(y,A,At,params,errx);
nmse_se1 = SE_DTurboCS(y,x,params);
nmse_se2 = SE_AMP(y,x,params);

% plot nmse
semilogy(1:iter,errx_out1,'r -x',1:iter,errx_out3,'k -*',1:iter,errx_out2,'g -->');
hold on;
semilogy(1:iter,nmse_se1,'b o',1:iter,nmse_se2,'m -<');
legend('Simulation, LET-Turbo-CS','Simulation, MMSE-Turbo-CS','Simulation, LET-AMP','State Evolution, LET-Turbo-CS','State Evolution, LET-AMP');
xlabel('Iteration');
ylabel('NMSE');
title('Comparisons of algorithms for sparse signal recovery');



