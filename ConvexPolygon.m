classdef ConvexPolygon
    % a 2D convex shape defined by points
    
    properties
        matrix_orig
        matrix
        % x1 y1
        % x2 y2
        % ...
        % xN yN
    end
        
    methods
        function obj = ConvexPolygon(mat)
            obj.matrix_orig = mat;
            obj.matrix = mat(grahamscan(mat),:);
        end

        function point = support(obj,direction)
            % input direction is a column
            % output point is a column
            perpdir = [-direction(2);direction(1)];
            dotprod=obj.matrix*direction;
            dotprod_perp = obj.matrix*perpdir;
            maxval = max(dotprod);
            paralelcandidates = find(dotprod==maxval);
            [~,argmax] = max(dotprod_perp(dotprod==maxval));
            point = obj.matrix(paralelcandidates(argmax),:).';
        end
        
        function mat = convexify(mat)
            mat = grahamscan(mat);
        end
        
    end
end
