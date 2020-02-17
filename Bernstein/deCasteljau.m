function b = deCasteljau(p,i,j,t)
    if i == 1
        b = p(j);
    else
        b = (1-t) .* deCasteljau(p,i-1,j-1,t) + t .* deCasteljau(p,i-1,j,t);
    end
end