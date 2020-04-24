function mon = BernsteinToMon(p, T)
    
    if size(p,1)==1 && size(p,2) ~= 1
        p = p';
    end
    
    N = size(p,1)-1;
    mon = zeros(size(p));
    for k = 0:N
        for i = 0:k
            mon(k+1,:) = mon(k+1,:) + nchoosek_mod(N,k)*nchoosek_mod(k,i)*(-1)^(k-i)*p(i+1,:);
        end
    end
    mon = flipud(mon);
    
    mon = mon ./ ( T.^(N:-1:0)');
end
