% compressive image recovery algroithm based on D-Turbo-CS
function [errout,x_hat]=TurboCS_img(y,A,At,params,errx)
% input signal: y
% measurement linear operators: A, At
%-------------parameters----------------
M=params.M;
N=params.N;
sigma=params.sigma;
K=params.K;
%-----------
delta=M/N;
iter=20;
epsilon=1e-6;
x_hat=zeros(N,1);
v = norm(y-A(x_hat))^2/M;
%-----------iterative algorithm------------------
for ii=1:iter
    r_hat = x_hat+1/delta*At(y-A(x_hat));
    c = (1/delta-1)*v+1/delta*sigma^2;
    % denoisers
    if K==3
        [F,F_div]=Kernel_lin_1(r_hat,c); 
    elseif K==2
        [F,F_div]=Kernel_lin_3(r_hat,c);
    elseif K==1
        [F,F_div]=BM3D_denoiser(r_hat,sqrt(c));
    end
    F_kernel=zeros(K,N);
    for jj=1:K
        F_kernel(jj,:)=F(jj,:)-F_div(jj)*reshape(r_hat,1,N);
    end
    %---sure minimization----
    C_ff=F_kernel*F_kernel'/N+epsilon*eye(K);
    C_fy=conj(F_kernel)*r_hat/N;
    C_cal=mldivide(C_ff,C_fy);
    % estimated signal
    x_hat=F_kernel.'*C_cal;
    v = norm(y-A(x_hat))^2/M;
    %vv  = 1/N*norm(x_hat-r_hat)^2-c;
    errout(ii) = errx(x_hat); % stats
end
%-----------output----------
P_ff=F*F'/N;
P_ff=P_ff+epsilon*eye(K);
P_fx=1/N*F*r_hat-c*F_div';
a=mldivide(P_ff,P_fx);
x_hat=F'*a;
end