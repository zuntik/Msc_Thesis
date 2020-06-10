close all;

p1 = [0 0; 0 1; 1 1; 1 0];
p2 = [2 0; 2 1; 3 1; 3 0];

p1 = [0 0; 0 1; 1 1; 1 0];
p2 = p1+[0.5 0.5];

% p1 = [4 11; 9 9; 4 5];
% p2 = [5 7; 12 7; 10 2; 7 3];

shape1=ConvexPolygon(p1);
shape2=ConvexPolygon(p2);

[intersection,simplex] = gjk2(shape1,shape2);
disp(intersection)

figure,hold on,axis equal
plot(polyshape(p1))
plot(polyshape(p2))

x = arrayfun(@(x)simplex{x}.P(1),1:length(simplex));
y = arrayfun(@(x)simplex{x}.P(2),1:length(simplex));
plot(polyshape(x,y))
scatter(0,0,'filled')

