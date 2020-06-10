function [intersection,simplex] = gjk2(shapeA,shapeB)
%function [dist, pt1, pt2, simplex] = gjk2(shapeA,shapeB)

    %[S,pts] = supportwrapper(shape1,shape2,[1 0].');
    % point is a struct
    d=[1;0];
    simplex={supportwrapper(shapeA,shapeB,d)};

    d = -d; 
    
    while true
        simplex{end+1} = supportwrapper(shapeA,shapeB,d) ; %#ok<AGROW>
        if(dot(simplex{end}.P,d)<=0)
            intersection=false;
            break
        end
        [intersection,simplex,d] = containsorigin(simplex);
        if intersection
            break;
        end
    end

end


function [intersection,simplex,d] = containsorigin(simplex)
    a = simplex{end};
    a0 = -a.P;
    if length(simplex) == 3
        b = simplex{end-1};
        c = simplex{end-2};
        ab = b.P - a.P;
        ac = c.P - a.P;
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
        ab = b.P - a.P;
        abPerp = tripleProduct(ab,a0,ab);
        d = abPerp;
        if ~any(d)
            intersection = true;
            return
        end
            
    end

    intersection=false;
end


function [point] = supportwrapper(shapeA,shapeB,D)
    % shapeA and shapeB are structs
    % D is a col
    point = struct();
    point.dir=D;
    ptA = shapeA.support(D);
    ptB = shapeB.support(-D);
    point.A=ptA;
    point.B=ptB;
    point.P=point.A-point.B;
    point.mag=sum(point.P.^2,'all');
end


function [vec] = tripleProduct(A,B,C)
    vec = B.*(dot(C,A))-A.*(dot(C,B));
end


function [dist,simplex] = mindist(simplex,shapeA,shapeB)
    
    d = closestpointto0(simplex{end}.P,simplex{end-1}.P);

    while true
        d = -d;
        c = supportwrapper(shapeA,shapeB,d);
        dc = dot(c.P,d);
        da = dot(d,a);
        if dc-da<10e-6
            dist = norm(d);
            return
        end
        p1 = closestpointto0(simplex{end}.P,c.P);
        p2 = closestpointto0(c.P,simplex{end-1},P);
        
        if p1(1).^2+p1(2).^2 < p1(1).^2+p2(2)^2
            simplex{end-1}= c;
            d = p1;
        else
            simplex{end} = c;
            d = p2;
        end
    end

    
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
