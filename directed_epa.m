function vec = directed_epa(support, simplex,d)
    % support is the support vector for the Minkowski differnece
    % simplex contains the origin and is of maximum dim (3)
    % d is direction to move
    
    extrabit = 10e-6;
    dPerp = [ -d(2), d(1) ];
    
    % edge 1
    for i = 1:4
        if i == 4
            disp('error, the origin is either on the edge or ouside the simplex');
        elseif segrayintersect([0;0],d,simplex{i}.vert,simplex{rem(i,3)+1}.vert)
            A = simplex{i}.vert;
            B = simplex{rem(i,3)+1}.vert;
            break
        end
    end
        
    tripleProduct = @(A,B,C) B.*(dot(C,A))-A.*(dot(C,B));
    while true
        e = B-A;
        oa = A;
        n = tripleProduct(e,oa,e);
        C = support(n);
        if abs(dot(C,dPerp)) < 0.000001
            vec = C + extrabit;
            return
        elseif all(C == B) || all(C == A)
            break
        elseif segrayintersect([0;0],d,A,C)
            B = C;
        elseif segrayintersect([0;0],d,B,C)
            A = C;
        else
            disp('error. at least one of these edges should intersect ray');
        end
    end
    
    % found the segment of the difference, now find the intersection
    [~,t,~] = segrayintersect([0;0],d,A,B);
    vec = d.*t;
    
    vec = vec + extrabit;
end


