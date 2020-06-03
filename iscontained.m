function relation = iscontained(polygon1,polygon2)
% returns true if polygon 2 is inside polygon 1

indexes = grahamscan([polygon1;polygon2]);
relation =  ~any(indexes(size(polygon1,1)+1:end));

end