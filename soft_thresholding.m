%--------------------------------------------------------------------------
function [F,F_div] = soft_thresholding(y, v)
alpha_2 = 2 * sqrt(v); % threshold set to 2*sqrt(v)

N = length(y);
y = reshape(y, 1, N);
%---------------------------------------------

index_1 = find( y < - alpha_2);
index_2 = find( y > alpha_2 );

f1 = zeros(1,N);
f1(index_1) = y(index_1) + alpha_2;
f1(index_2) = y(index_2) - alpha_2;
F = f1;

eta=randn(1,N);
epsilon=max(y(:))/1000+eps;
y2=y+epsilon*eta;

index_1 = find( y2 < - alpha_2);
index_2 = find( y2 > alpha_2 );

f2 = zeros(1,N);
f2(index_1) = y2(index_1) + alpha_2;
f2(index_2) = y2(index_2) - alpha_2;
F2 = f2;

F_div=eta*((F2-F)/epsilon)'/N;


