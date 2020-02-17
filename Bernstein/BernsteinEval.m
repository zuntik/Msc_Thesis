function [points] = BernsteinEval(p,T,times)
    [num_points,dim] = size(p);
    points = zeros(dim,length(times));    
    for t = 1:length(times)
        for i = 1:dim
            points(i,t) = deCasteljau(p(:,i), num_points, num_points, times(t)/T);
        end
    end
end