s1 = [ 4 11; 9 9;  4 5 ];
s2 = [ 5 7;  12 7; 10 2; 7 3 ];

shape1 = ConvexPolygon(s1);
shape2 = ConvexPolygon(s2);

d = [0; 1];

%[dist, pt1, pt2, simplex] = gjk2(shape1,shape2);
[intersect,simplex] = gjk2(shape1,shape2);

disp(intersect)

if intersect && false
    disp(directed_epa(@(d)shape1.support(d)-shape2.support(-d),simplex,d))
end

box1 = [ 0 0; 0 1 ; 1 0; 1 1 ];
box2 = box1 + [ 20 0 ];
shape3 = ConvexPolygon(box1);
shape4 = ConvexPolygon(box2);

[intersect,simplex] = gjk2(shape3,shape4);
disp(intersect)
