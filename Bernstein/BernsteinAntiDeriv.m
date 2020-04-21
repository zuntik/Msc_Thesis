function [new_p] = BernsteinAntiDeriv(p,T,p0)
    [ n, dim ] = size(p);
    if n==1 && dim ~= 1
        p = p';
        [ n, dim ] = size(p);
    end
    
    new_p = zeros(n+1,dim);
    for i = 1:size(p,1)
        new_p(i+1,:) = sum(p(1:i,:),1);
    end
    new_p = p0 + new_p.*(T/size(p,1));
end