function [intersect, tr, ts] = segrayintersect(rayOrigin, rayDirection, point1,  point2)
    intersect = false;
    v1 = rayOrigin - point1;
    v2 = point2 - point1;
    v3 = [-rayDirection(2) ; rayDirection(1) ];
    dotprod = dot(v2,v3);
    if abs(dotprod) < 0.000001
        reuturn
    end

    crossprod = @(a,b)  (a(1)*b(2)) - (a(2)*b(1));
    tr = crossprod(v2,v1) / dotprod;
    ts = dot(v1,v3) / dotprod;
    if tr >= 0.0 && ts >= 0.0 && ts <= 1.0
        intersect = true;
    end

end