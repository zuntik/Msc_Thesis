classdef ConvexPolygon
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
            point = obj.matrix(argmax,:).';
        end
    end
end
