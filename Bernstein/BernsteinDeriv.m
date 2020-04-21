function [new_p] = BernsteinDeriv(p,T)
    if size(p,1)==1 && size(p,2) ~= 1
        p = p';
    end
    new_p = (p(2:end,:)-p(1:end-1,:)).*((size(p,1)-1)/T);
end