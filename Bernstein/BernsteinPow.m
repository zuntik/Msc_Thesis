function new_p = BernsteinPow(p, y)
    if y == 0
        new_p = ones(1,size(p,2));
        return
    end
    temp_p = BernsteinPow(p,floor(y/2));
    if rem(y,2) == 0
        new_p = BernsteinMul(temp_p,temp_p);
    else
        new_p = BernsteinMul(p,BernsteinMul(temp_p,temp_p));
    end
end
