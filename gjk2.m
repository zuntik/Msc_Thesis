function [intersection,simplex] = gjk2(shapeA,shapeB)
%function [dist, pt1, pt2, simplex] = gjk2(shapeA,shapeB)

    %[S,pts] = supportwrapper(shape1,shape2,[1 0].');
    % point is a struct
    point = supportwrapper(shapeA,shapeB,[1 0].');
    
    %simplex = {struct('vert',S,'pts',pts)};
    simplex = {struct('vert',point)};
    d = [-1 0].'; 

    while true
        simplex{end+1} = struct('vert',supportwrapper(shapeA,shapeB,d));
        %[A,pts] = supportwrapper(shapeA,shapeB,D);
        if dot(simplex{end}.vert,d) < 0
            intersection = false;
            break
        end
        %simplex{end+1} = struct('vert',A,'pts',pts);
        [intersection,simplex,d] = containsorigin(simplex);
        if intersection
            break;
        end
    end

end


function [intersection,simplex,d] = containsorigin(simplex)
    a = simplex{end};
    a0 = -a.vert;
    if length(simplex) == 3
        b = simplex{end-1};
        c = simplex{end-2};
        ab = b.vert - a.vert;
        ac = c.vert - a.vert;
        abPerp = tripleProduct(ac,ab,ab);
        acPerp = tripleProduct(ab,ac,ac);
        if dot(abPerp,a0)> 0
            simplex(end-2) = [];
            d = abPerp;
        elseif dot(acPerp,a0) > 0
            simplex(end-1) = [];
            d = acPerp;
        else
            intersection=true;
            d = [];
            return
        end
    else
        b = simplex{1};
        ab = b.vert - a.vert;
        abPerp = tripleProduct(ab,a0,ab);
        d = abPerp;
    end

    intersection=false;
end


function [diff,pts] = supportwrapper(shapeA,shapeB,D)
    % shapeA and shapeB are structs
    % D is a col
    % A is a col
    % pts are 2 cols side by side
    ptA = shapeA.support(D);
    ptB = shapeB.support(-D);
    diff = ptA-ptB;
    pts = [ ptA,ptB ];
end


function [vec] = tripleProduct(A,B,C)
    vec = B.*(dot(C,A))-A.*(dot(C,B));
end


function [closestPt, t] = dist2line(A, B)
%DIST2LINE Finds the closest point on a line segment to the origin.
%
%   INPUTS
%   A: N dimensional point on the line segment AB.
%   B: N dimensional point on the line segment AB.
%
%   RETURNS
%   closestPt: N-D dimensional point on line segment AB which is closest to
%       the origin.
%   t: Weight used to solve for the point on the Bernstein polynomial,
%       where f(t) = (1-t)*A + t*B, t = [0, 1]

    if all(abs(A - B) < 1000*eps)
        t = 0;
        closestPt = A;
        return
    end

    v = B - A;
    u = A;
    t = -dot(v,u)/dot(v,v);

    % If t is out of the range [0, 1] we are closest to A or B
    if t > 1
        t = 1;
        closestPt = B;
    elseif t < 0
        t = 0;
        closestPt = A;
    else
        closestPt = (1-t)*A + t*B;
    end
end


function [intersects, tr,ts] = segrayintersect(rayOrigin, rayDirection, point1, point2)
% https://stackoverflow.com/a/32146853

    v1 = rayOrigin - point1;
    v2 = point2 - point1;
    v3 = [-rayDirection(2), rayDirection(1)].';

    dot_prod = dot(v2,v3);
    if abs(dot_prod) < 0.000001
        intersects = false;
        return
    end
    
    cross_prod = @(v1,v2) (v1(1)*v2(2)) - (v1(2)*v2(1));
    tr = cross_prod(v2,v1) /dot_prod;
    ts = dot(v1,v3)/dot_prod;
    if (tr >= 0.0 && (ts >= 0.0 && ts <= 1.0))
        intersects = true;
    else
        intersects = false;
    end
end
