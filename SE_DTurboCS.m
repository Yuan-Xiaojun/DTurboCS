% state evolution of D-Turbo-CS algorithm for random partial orthogonal linear
% mearuement
function nmse = SE_DTurboCS(y,x,params)
N = params.N;
M = params.M;
K = params.K;
sigma = params.sigma;
iter = params.iter;
type = params.type;
v = norm(y)^2/M-sigma^2;
epsilon=1e-6;
for ii=1:iter
    % for partial orthogonal linear operator
    if type=="orthogonal"
        c = (N/M-1)*v+N/M*sigma^2;
    else
        c = (N/M)*v+N/M*sigma^2;
    end
    sumtheta=0;
    sumtheta2=0;
    for kk=1:10
        r=x+sqrt(c)*randn(N,1);
        
        if K==3
            [F,F_div]=Kernel_lin_1(r,c);
        elseif K==2
            [F,F_div]=Kernel_lin_3(r,c);
        elseif K==1
            [F,F_div]=BM3D_denoiser(r,sqrt(c));
        end
        F_kernel=F-F_div'*r';
        C_ff=F_kernel*F_kernel'/N+epsilon*eye(K);
        C_fy=conj(F_kernel)*r/N;
        C_cal=mldivide(C_ff,C_fy);
        x_state=F_kernel.'*C_cal;
        
        %
        P_ff=F*F'/N;
        P_ff=P_ff+epsilon*eye(K);
        P_fx=1/N*F*r-c*F_div';
        a=mldivide(P_ff,P_fx);
        x_hat=F'*a;
        
        theta2=1/N*norm(x_state-x)^2;
        theta22=1/N*norm(x_hat-x)^2;
        sumtheta=sumtheta+theta2;
        sumtheta2=sumtheta2+theta22;
    end
    
    v=sumtheta/10;
    v2=sumtheta2/10;
    nmse(ii) = v2;
end
end