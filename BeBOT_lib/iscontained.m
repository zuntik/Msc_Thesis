function relation = iscontained(polygon1,polygon2)
% returns true if polygon 2 is inside polygon 1
polygon1 = polygon1(:,1:2);
polygon2 = polygon2(:,1:2);
indexes = grahamscan([polygon1;polygon2]);
relation =  ~any(indexes(size(polygon1,1)+1:end));

end