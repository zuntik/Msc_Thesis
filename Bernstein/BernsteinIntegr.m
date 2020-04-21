function [val] = BernsteinIntegr(p,T)
    if size(p,1)==1 && size(p,2) ~= 1
        p = p';
    end
    val = T / size(p,1) * sum(p,1);
end