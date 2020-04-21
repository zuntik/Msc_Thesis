function d = GJK2(shape1,shape2)
% GJK2
    dim = shape1.dim;
    if shape2.dim ~= dim
        disp('error! shape1 and shape2 have different dimensions');
    end

    function difference = support(s1, s2, d)
        p1 = s1.getFarthestPointInDirection(d);
        p2 = s2.getFarthestPointInDirection(-d);
        difference = p1 - p2;
    end

    simplex = Simplex();

    d = zeros(dim,1);
    d(1) = 1;

    while True
        simplex.add(support(s1,s2,d));
        if simplex.getLast().dot(d) <=0 
            return False
        end
    end

end