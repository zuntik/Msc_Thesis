function [points] = BernsteinEval(p,T,times)
% a function to evaluate a bezier curve with control points given by p
%   p contains the control points in the following example form:
%     x1 y1 
%     x2 y2
%     x3 y3
%   T is the time length of the curve
%   times are the points to evaluate

    [num_points,dim] = size(p);
    
    if num_points==1 && dim ~= 1
        p = p';
        dim = dim+num_points;
        num_points= dim - num_points;
        dim = dim- num_points;
    end
    
    points = zeros(dim,length(times));    
    for t = 1:length(times)
        for i = 1:dim
            points(i,t) = deCasteljau(p(:,i), num_points, num_points, times(t)/T);
        end
    end
end