function [mat] = BernsteinEvalMat(n,T,times)
% a function to evaluate a bezier curve with control points given by p
%   n is the order (n+1 control points)
%   T is the time length of the curve
%   times are the points to evaluate

    mat = zeros(length(times),n+1);
    for i = 1:length(times)
        for j = 0:n
            mat(i,j+1) = BernsteinBasis(j,n,times(i)/T);
        end
    end
end