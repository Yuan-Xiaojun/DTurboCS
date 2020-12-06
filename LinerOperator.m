%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear operators of different types
% the l2-norm of the rows of A are normalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [A, At] = LinerOperator(n,m,type)
    perm=randperm(n);
    indexs=perm(1:m);
    sign1=2*(rand(n,1)>0.5)-1;
    sign2=2*(rand(n,1)>0.5)-1;
    if(type=="highrandsecdct")
        A=@(z)subsref(idct(sign1.*dct(sign2.*z(:))),struct('type','()','subs',{{indexs}}));%F'S1FS2
        At=@(z)reshape(sign2.*idct(sign1.*dct(put_vector(n,indexs,z))),[n,1]);
    end
    if(type=="highfftrandsec")
        A=@(z)subsref(sqrt(n)*ifft(1/sqrt(n)*sign1.*fft(sign2.*z(:))),struct('type','()','subs',{{indexs}}));%F'S1FS2
        At=@(z)reshape(sqrt(n)*sign2.*ifft(1/sqrt(n)*sign1.*fft(put_vector(n,indexs,z))),[n,1]);
    end
    if(type=="randsecdctsign")
        A=@(z)subsref(dct(sign2.*z(:)),struct('type','()','subs',{{indexs}}));%F'S1FS2
        At=@(z)reshape(sign2.*idct(put_vector(n,indexs,z)),[n,1]);
    end
    
    if(type=="randsecdct")
        A=@(z) subsref(dct(z(:)),struct('type','()','subs',{{indexs}}));
        At=@(z) reshape(idct(put_vector(n,indexs,z)),[n,1]);
    end
    if(type=="gaurand")
        GA=randn(m,n)/sqrt(n);
        A=@(z) GA*z(:);
        At=@(z) reshape(GA'*z,[n,1]);
    end
    if(type=="mc")
         A=@(z) subsref(z(:),struct('type','()','subs',{{indexs}}));
         At=@(z) reshape(put_vector(n,indexs,z),[n,1]);
    end
end