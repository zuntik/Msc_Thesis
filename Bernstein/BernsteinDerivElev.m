function [new_p] = BernsteinDerivElev(p,T)
    if size(p,1)==1 && size(p,2) ~= 1
        p = p';
    end
    new_p = BernsteinDeriv(BernsteinDegrElev(p,size(p,1)),T);
    %new_p = BernsteinDegrElev(BernsteinDeriv(p,T),size(p,1));
end