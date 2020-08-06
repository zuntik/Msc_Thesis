function val = logbarrierfunc(delta,z)
    z = my_sigma(z);
    k = 2;
    if z>delta
        val = - log(z);
    else
        val = ((k-1)/k) * (((z-k*delta)/((k-1)*delta))^k-1) - log(delta);
        %val = 0.5*(((z-2)).^2 -1)-log(delta);
    end
end

function val = my_sigma(z) %#ok<DEFNU>
    if z >= 0
        val = tanh(z);
    else
        val = z;
    end
end