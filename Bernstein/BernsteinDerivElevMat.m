function mat = BernsteinDerivElevMat(n,T)
	% Returns matrix used for derivation and elevation
	% INPUT
	% n: order (= nº control points  - 1)
	% T: time
	% OUTPUT
    % mat
    
    mat = BernsteinDerivMat(n+1,T)*BernsteinDegrElevMat(n,n+1);
    %mat = BernsteinDegrElevMat(n-1,n)*BernsteinDerivMat(n,T);
    
end