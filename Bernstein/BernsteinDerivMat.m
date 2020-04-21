function [mat] = BernsteinDerivMat(n,T)
%BERNSTEINDERIVMAT return matrix to derivate control points
%   n - order
%   T - time

mat = n/T .*([zeros(n,1),eye(n)]-[eye(n),zeros(n,1)]);

end

