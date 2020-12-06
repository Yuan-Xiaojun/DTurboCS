%--------------------------------------------------------------------------
% third kind of piece-wise linear function. See
% "Near Optimal compressed sensing without priors: parametric SURE approximate
% message passing"
function [F,F_div] = Kernel_lin_3(y,v)

T=2*sqrt(v);

n=length(y);
f1=zeros(1,n);
f2=zeros(1,n);

f1=y;
f2=y.*exp(-y.^2/(2*T^2));

F(1,:) = f1;
F(2,:) = f2;

F_div(1)=1;
F_div(2)=1/n*sum(exp(-y.^2/(2*T^2)).*(1-y.^2/T^2));
end