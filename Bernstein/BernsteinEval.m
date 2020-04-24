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
    end
    [num_points,dim] = size(p);
    if dim==1
        points = arrayfun(@(t)BernsteinEvalMat(size(p,1)-1,T,t)*p,times);
    else
        points = BernsteinEvalMat(size(p,1)-1,T,times(:))*p;
    end
    
end