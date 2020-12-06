% state evolution of D-Turbo-CS algorithm for random partial orthogonal linear
% mearuement
function nmse = SE_AMP(y,x,params)
N = params.N;
M = params.M;
K = params.K;
sigma = params.sigma;
iter = params.iter;
c = norm(y)^2/M;
epsilon=1e-6;
for ii=1:iter
    % for partial orthogonal linear operator
    if ii~=1
        c = v*(N/M)+(N/M)*sigma^2;
    end
    sumtheta=0;
    for kk=1:10
        r=x+sqrt(c)*randn(N,1);
        if K==3
            [F,F_div]=Kernel_lin_1(r,c);
        elseif K==2
            [F,F_div]=Kernel_lin_3(r,c);
        elseif K==1
            [F,F_div]=BM3D_denoiser(r,sqrt(c));
        end
        P_ff=F*F'/N+epsilon*eye(1);
        P_fx=1/N*F*r-c*F_div';
        a=mldivide(P_ff,P_fx);

        x_hat=F'*a;
        theta2=1/N*norm(x_hat-x)^2;
        sumtheta=sumtheta+theta2;
    end
    v=sumtheta/10;
    nmse(ii) = v;
end
end