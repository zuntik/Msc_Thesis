function [new_p] = BernsteinAntiDeriv(p,T,p0)
    new_p = zeros(size(p,1)+1,size(p,2));
    for i = 1:size(p,1)
        new_p(i+1,:) = sum(p(1:i,:),1);
    end
    new_p = p0 + new_p.*(T/size(p,1));
end