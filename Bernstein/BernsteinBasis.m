function val = BernsteinBasis(k,n,tau)
% n is order
    val = nchoosek_mod(n,k)*(1-tau)^(n-k)*tau^k;
end