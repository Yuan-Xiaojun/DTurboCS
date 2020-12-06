% mmse estimator for Gaussian-Bernoulli distributed sparse signal
% lambda: sparsity
function [xpost,vpost,xext,vext] = mmse_denoiser(xpri,vpri,lambda)
    EXP_MAX = 50;
    EXP_MIN = -50;
    
    exponent =  - (abs(xpri).^2 * 1/lambda / ( vpri * (vpri+1/lambda)))/2;
    exponent(exponent > EXP_MAX) = EXP_MAX;
    exponent(exponent < EXP_MIN) = EXP_MIN;
            
    den = 1 + sqrt(( vpri + 1/lambda ) / vpri) * (1-lambda) / lambda * exp(exponent);
    C = 1./den;
            
    % compute a posteriori
    xpost = ( 1/lambda ) / ( 1/lambda + vpri ) * xpri .* C;
    VAR = 1/lambda * vpri / ( vpri + 1/lambda ) * C + abs( 1/lambda / ( 1/lambda + vpri ) * xpri ).^2 .* C .* (1-C);
    vpost = mean(VAR);
            
    % update extrinsic
    vext = 1/(1/vpost-1/vpri);
    xext = ( vext / vpost ) * xpost - (vext / vpri) * xpri;
end