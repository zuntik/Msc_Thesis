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
        elseif segrayintersect([0;0],d,simplex{i}.P,simplex{rem(i,3)+1}.P)
            A = simplex{i}.P;
            B = simplex{rem(i,3)+1}.P;
            break
        end
    end
    
    % if oritentation > 0 clockwise, elseif 0 colinear else anticlockwise
    orientation = @(p1,p2,p3) (p2(2)-p1(2))*(p3(1)-p2(1))-(p2(1)-p1(1))*(p3(2)-p2(2));
    if orientation(A,B,[0 0]) <= 0
        % swap
        A = A+B;
        B = A-B;
        A = A-B;
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


