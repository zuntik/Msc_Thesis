function [val] = BernsteinIntegr(p,T)
    val = T / size(p,1) * sum(p,1);
end