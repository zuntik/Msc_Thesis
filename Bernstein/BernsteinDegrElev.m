function new_cp = BernsteinDegrElev(cp, m)
	% Performs degree elevation of curve defined by cp to M
	% INPUT
	% cp: matrix N x dim
	% M: new desired order
	% OUTPUT
	% new_cp: new control points for new order

	[num_points,dim] = size(cp);
	
	if num_points==1 && dim ~= 1
		cp = cp';
		dim = dim+num_points;
		num_points= dim - num_points;
		dim = dim- num_points;
    end
    
    n = num_points-1;
    
    if n==m
        return
    end
    assert(m>n,'new degree not greater than current');
    
	new_cp = zeros(m+1, dim);
	
    for i = 0:m
		for j = [ max(0, i-m+n) : min(n,i) ]
			%new_cp(i+1,:) = new_cp(i+1,:) + (nchoosek_mod(n,j)*nchoosek_mod(m-n,i-j)/nchoosek_mod(m,i)).*cp(j+1,:);
			new_cp(i+1,:) = new_cp(i+1,:) + nchoosek_mod(i,j)*nchoosek_mod(m-i,n-j).*cp(j+1,:);
		end
    end
    new_cp = new_cp./nchoosek_mod(m,n);

end