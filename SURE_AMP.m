function [nmse,x_hat]=SURE_AMP(y,A,At,params,errx)
% input signal: y
% measurement matrix: A, At
%-------------parameters----------------
M=params.M;
N=params.N;
K=params.K;
delta=M/N;
iter=params.iter;
x_hat=zeros(N,1);
xr = params.xr; % real value of x
type = params.type;
sigma=params.sigma;
epsilon=1e-6;

z = y*sqrt(N/M);
c = norm(y)^2/M+N/M*sigma^2;
%-----------iterative algorithm------------------
for ii=1:iter
    r_hat = (x_hat+sqrt(N/M)*At(z));
    if type=="orthogonal"
        if ii~=1
            c=1/log(2)*median(abs(r_hat))^2;
        end
    else
        c = norm(z)^2/M;
    end
    % Kernel function to promote sparsity
    if K==3
        [F,F_div]=Kernel_lin_1(r_hat,c); 
    elseif K==2
        [F,F_div]=Kernel_lin_3(r_hat,c);
    elseif K==1
        [F,F_div]=BM3D_denoiser(r_hat,sqrt(c));
    end

    P_ff=F*F'/N+epsilon*eye(1);
    P_fx=1/N*F*r_hat-c*F_div';
    a=mldivide(P_ff,P_fx);
       
    x_hat=F'*a;
    v=F_div*a;
    z=sqrt(N/M)*(y-A(x_hat))+1/delta*v*z;
    
    nmse(ii) = errx(x_hat);
end
end