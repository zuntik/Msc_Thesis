function [dist, pt1, pt2, simplex] = gjk2(shape1,shape2)

    [S,pts] = supportwrapper(shape1,shape2,[1 0].');

    simplex = {struct('vert',S,'pts',pts)};
    D = -S; 

    while True
        [A,pts] = supportwrapper(shapeA,shapeB,D);
        if dot(A,D) < 0
            intersection = False;
            break
        end
        simplex{end+1} = struct('vert',A,'pts',pts);
        [intersection,simplex,D] = doSimplex(simplex);
        if intersection
            break
        end
    end



end


function simplex = doSimplex(simplex)
    switch length(simplex)
        case 2
            simplex = linesimplex(simplex);
        case 3
            simplex = triagsimplex(simplex);
    end
end

function simplex = linesimplex(simplex)
    if 
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


classdef ConvexFigure 
    methods(Abstract)
        point = support(direction);
    end
end


classdef ConvexPolygon < ConvexFigure
    % a 2D convex shape defined by points
    
    properties
        matrix
        % x1 y1
        % x2 y2
        % ...
        % xN yN
    end
        
    methods
        function obj = ConvexPolygon(mat)
            obj.matrix = mat;
        end
            
        function point = support(obj,direction)
            % input direction is a column
            % output point is a column
            [~,argmax] = max(obj.matrix * direction);
            point = obj(:,argmax).';
        end
    end
end