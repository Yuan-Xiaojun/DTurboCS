% do xs[index] = x, where xs is a vector of size n
function xs=put_vector(n,index,x)
xs=zeros(n,1);
xs(index)=x;
end