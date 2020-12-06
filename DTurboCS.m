function [nmse,x_hat]=DTurboCS(y,A,At,params,errx)
% input signal: y
% measurement matrix: A, At
%-------------parameters----------------
M=params.M;
N=params.N;
sigma=params.sigma;
K=params.K;
delta=M/N;
iter=params.iter;
x_hat=zeros(N,1);
xr = params.xr; % real value of x
type = params.type;
epsilon=1e-6;
%-----------iterative algorithm------------------
for ii=1:iter
    v = norm(y-A(x_hat))^2/M-sigma^2;
    r_hat = x_hat+1/delta*At(y-A(x_hat));
    if type=="orthogonal"
        c = (1/delta-1)*v+1/delta*sigma^2;
    else
        c = (1/delta)*v+1/delta*sigma^2;
    end
    % Kernel function to promote sparsity
    if K==3
        [F,F_div]=Kernel_lin_1(r_hat,c); 
    elseif K==2
        [F,F_div]=Kernel_lin_3(r_hat,c);
    elseif K==1
        [F,F_div]=BM3D_denoiser(r_hat,sqrt(c));
    end
    
    F_kernel=F-F_div'*r_hat';
    %---estimate of x by sure minimization----
    C_ff=F_kernel*F_kernel'/N+epsilon*eye(K);
    C_fy=conj(F_kernel)*r_hat/N;
    C_cal=mldivide(C_ff,C_fy);
    % estimated signal
    x_hat=F_kernel.'*C_cal;
    %
    P_ff=F*F'/N;
    P_ff=P_ff+epsilon*eye(K);
    P_fx=1/N*F*r_hat-c*F_div';
    a=mldivide(P_ff,P_fx);
    x_hat2=F'*a;

    nmse(ii) = errx(x_hat2);
end
%-----------sure let output step----------
P_ff=F*F'/N;
P_ff=P_ff+epsilon*eye(K);
P_fx=1/N*F*r_hat-c*F_div';
a=mldivide(P_ff,P_fx);
x_hat=F'*a;
end